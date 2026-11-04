% FUSELAGE_BASELINE_STUDY  Baseline laminate family comparison for a fuselage skin panel.
% Illustrative, preliminary-design-level loads and strengths -- NOT certified values.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));

mat = struct('E1', 161e9, 'E2', 11.4e9, 'G12', 5.17e9, 'nu12', 0.32, ...
             'density', 1580.0, 'thickness', 0.184e-3);
S   = struct('Xt', 2800e6, 'Xc', 1500e6, 'Yt', 60e6, 'Yc', 200e6, 'S12', 90e6);

% Representative membrane load case (thin-wall pressurization + illustrative increments)
dp = 58000;  R = 2.0;                    % cabin differential pressure (Pa), skin radius (m)
Nx  = dp*R/2 + 40000;                    % N/m: end-cap equilibrium + representative bending
Ny  = dp*R;                              % N/m: hoop resultant
Nxy = 20000;                             % N/m: representative shear / torsion
N = [Nx; Ny; Nxy];  M = [0; 0; 0];       % membrane-only (no bending moments)

res = run_baseline_study(mat, S, N, M);
fprintf('%-34s %5s %7s %9s %8s %8s %8s %7s %7s %7s  %s\n', 'Laminate', 'Plies', ...
        't(mm)', 'm(kg/m2)', 'A11', 'A22', 'A66', 'MaxSt', 'TsaiH', 'TsaiW', 'Governing');
for i = 1:numel(res)
    r = res(i);
    fprintf('%-34s %5d %7.3f %9.3f %8.2f %8.2f %8.2f %7.3f %7.3f %7.3f  %s\n', ...
        r.name, r.plies, r.thickness_mm, r.areal_mass, r.A11, r.A22, r.A66, ...
        r.fi_max_stress, r.fi_tsai_hill, r.fi_tsai_wu, r.governing);
end
