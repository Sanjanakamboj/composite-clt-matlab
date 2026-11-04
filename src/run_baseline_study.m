function results = run_baseline_study(mat, S, N, M)
%RUN_BASELINE_STUDY Compare the four baseline symmetric laminate families under a load.
%   Returns a struct array with name, angles (full stack), plies, thickness_mm,
%   areal_mass, A11/A22/A66 (MN/m), worst-case failure indices over all ply
%   surfaces, and the governing ply (Max Stress).
    names = {'quasi_isotropic_[0/45/-45/90]_s', 'zero_biased_[0/0/45/-45]_s', ...
             'ninety_biased_[90/90/45/-45]_s',  'pm45_biased_[45/-45/45/-45]_s'};
    halves = {[0 45 -45 90], [0 0 45 -45], [90 90 45 -45], [45 -45 45 -45]};
    results = struct('name', {}, 'angles', {}, 'plies', {}, 'thickness_mm', {}, ...
                     'areal_mass', {}, 'A11', {}, 'A22', {}, 'A66', {}, ...
                     'fi_max_stress', {}, 'fi_tsai_hill', {}, 'fi_tsai_wu', {}, ...
                     'governing', {}, 'pass', {});
    for i = 1:numel(names)
        angles = [halves{i}, fliplr(halves{i})];          % "_s" symmetric stack
        [ABD, z, am] = laminate_abd(mat, angles);
        [eps0, kappa] = solve_deformation(ABD, N, M);
        st = recover_ply_states(mat, angles, z, eps0, kappa);
        ms = 0; th = -inf; tw = -inf; gov = 1;
        for j = 1:numel(st)
            [fi, ~] = failure_indices(st(j).material_stress, S);
            if fi.max_stress > ms * (1 + 1e-9), ms = fi.max_stress; gov = j; end   % ties -> first listed
            th = max(th, fi.tsai_hill);
            tw = max(tw, fi.tsai_wu);
        end
        [~, gmode] = failure_indices(st(gov).material_stress, S);
        results(i).name         = names{i};
        results(i).angles       = angles;
        results(i).plies        = numel(angles);
        results(i).thickness_mm = numel(angles) * mat.thickness * 1e3;
        results(i).areal_mass   = am;
        results(i).A11 = ABD(1,1)/1e6;  results(i).A22 = ABD(2,2)/1e6;  results(i).A66 = ABD(3,3)/1e6;
        results(i).fi_max_stress = ms;  results(i).fi_tsai_hill = th;  results(i).fi_tsai_wu = tw;
        results(i).governing = sprintf('%+d deg (%s, %s)', st(gov).theta_deg, st(gov).surface, gmode);
        results(i).pass = max([ms th tw]) < 1;
    end
end
