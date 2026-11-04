function run_tests()
% RUN_TESTS  Verification of the MATLAB port against the Python original
% (composite-fuselage-laminate-analysis) and against closed-form checks.
% Works in MATLAB and GNU Octave.  Run from anywhere:  run_tests
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));

mat = struct('E1', 161e9, 'E2', 11.4e9, 'G12', 5.17e9, 'nu12', 0.32, ...
             'density', 1580.0, 'thickness', 0.184e-3);
S   = struct('Xt', 2800e6, 'Xc', 1500e6, 'Yt', 60e6, 'Yc', 200e6, 'S12', 90e6);
nfail = 0;

function check(name, ok)
    if ok, fprintf('  PASS  %s\n', name); else, fprintf('  FAIL  %s\n', name); end
end
function r = close_to(a, b, rtol)
    r = all(abs(a(:) - b(:)) <= rtol * max(abs(b(:)), eps) + 0);
end

% ---- Lamina stiffness (reference values from the Python implementation) ----
Q = reduced_stiffness(mat);
Qref = [162175885976.73166 3674643677.286442 0; 3674643677.286442 11483261491.52013 0; 0 0 5170000000];
ok = close_to(Q, Qref, 1e-12);  check('Q matches Python reference', ok); nfail = nfail + ~ok;

Qb = transformed_reduced_stiffness(mat, 30);
Qbref = [97197131084.11401 30980242448.60124 48390805042.78803;
         30980242448.60124 21850818841.508217 16861015440.783009;
         48390805042.78803 16861015440.783009 32475598771.314796];
ok = close_to(Qb, Qbref, 1e-12);  check('Qbar(30 deg) matches Python reference', ok); nfail = nfail + ~ok;

% ---- Closed-form identities ----
ok = close_to(transformed_reduced_stiffness(mat, 0), Q, 1e-12);
check('Qbar(0 deg) == Q', ok); nfail = nfail + ~ok;
Qp = transformed_reduced_stiffness(mat, 45); Qm = transformed_reduced_stiffness(mat, -45);
ok = close_to(Qp([1 2 5 9]), Qm([1 2 5 9]), 1e-12) && close_to(Qp(1,3), -Qm(1,3), 1e-12);
check('+45/-45: even terms equal, shear-coupling terms flip sign', ok); nfail = nfail + ~ok;

% ---- Strain/stress transform round trip (Reuter matrix) ----
e = [1e-3; -2e-3; 5e-4];  th = 37;
em = transform_strain_to_material(e, th);
Tinv = stress_transform_matrix(-th);
eback = diag([1 1 2]) * Tinv * diag([1 1 0.5]) * em;
ok = close_to(eback, e, 1e-12);  check('strain global->material->global round trip', ok); nfail = nfail + ~ok;

% ---- Laminate: symmetric => B = 0 ----
[ABD, ~, ~] = laminate_abd(mat, [0 45 -45 90 90 -45 45 0]);
ok = max(abs(ABD(1:3,4:6)), [], 'all') < 1e-6 * max(abs(ABD(:)));
check('symmetric laminate has B = 0', ok); nfail = nfail + ~ok;

% ---- Baseline study vs Python reference ----
N = [98000; 116000; 20000];  M = [0; 0; 0];
res = run_baseline_study(mat, S, N, M);
refA11 = [101017238.276 156472124.086 45562352.4652 74221344.0148] / 1e6;
refA22 = [101017238.276 45562352.4652 156472124.086 74221344.0148] / 1e6;
refA66 = [34406134.2609 34406134.2609 34406134.2609 61202028.5218] / 1e6;
refMS  = [0.24073179715312484 0.4784023531430562 0.3857806954910518 0.224170685867536];
refTH  = [0.05898555898331164 0.22998294942693606 0.15027093863763852 0.056152286039710196];
refTW  = [0.15535374617236006 0.3924022701573669 0.2905005911255662 0.13970858943177653];
ok = close_to([res.A11], refA11, 1e-9) && close_to([res.A22], refA22, 1e-9) && close_to([res.A66], refA66, 1e-9);
check('A11/A22/A66 for all 4 laminate families match Python', ok); nfail = nfail + ~ok;
ok = close_to([res.fi_max_stress], refMS, 1e-9);  check('Max-Stress FI matches Python', ok); nfail = nfail + ~ok;
ok = close_to([res.fi_tsai_hill], refTH, 1e-9);   check('Tsai-Hill FI matches Python', ok); nfail = nfail + ~ok;
ok = close_to([res.fi_tsai_wu], refTW, 1e-9);     check('Tsai-Wu FI matches Python', ok); nfail = nfail + ~ok;
ok = all([res.pass]);  check('all baseline laminates pass at this load', ok); nfail = nfail + ~ok;

% ---- Ill-conditioned ABD is rejected ----
threw = false;
try, solve_deformation(zeros(6), N, M); catch, threw = true; end
check('singular ABD raises an error', threw); nfail = nfail + ~threw;

if nfail == 0, fprintf('\nAll tests passed.\n'); else, fprintf('\n%d test(s) FAILED.\n', nfail); end
end

function check(name, ok)
    if ok, fprintf('  PASS  %s\n', name); else, fprintf('  FAIL  %s\n', name); end
end
function r = close_to(a, b, rtol)
    r = all(abs(a(:) - b(:)) <= rtol * max(abs(b(:)), eps));
end
