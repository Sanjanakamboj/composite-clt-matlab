function states = recover_ply_states(mat, angles_deg, z, eps0, kappa)
%RECOVER_PLY_STATES Strain/stress at the bottom and top surface of every ply.
%   Returns a struct array (2 entries per ply, bottom then top) with fields:
%   ply_index (1-based), theta_deg, z, surface, global_strain, material_strain,
%   global_stress, material_stress.
    Q = reduced_stiffness(mat);
    states = struct('ply_index', {}, 'theta_deg', {}, 'z', {}, 'surface', {}, ...
                    'global_strain', {}, 'material_strain', {}, ...
                    'global_stress', {}, 'material_stress', {});
    idx = 0;
    for k = 1:numel(angles_deg)
        Qbar = transformed_reduced_stiffness(mat, angles_deg(k));
        for s = 1:2
            if s == 1, zz = z(k);   surf = 'bottom';
            else,      zz = z(k+1); surf = 'top';
            end
            eg = eps0(:) + zz * kappa(:);              % linear through-thickness strain
            em = transform_strain_to_material(eg, angles_deg(k));
            idx = idx + 1;
            states(idx).ply_index       = k;
            states(idx).theta_deg       = angles_deg(k);
            states(idx).z               = zz;
            states(idx).surface         = surf;
            states(idx).global_strain   = eg;
            states(idx).material_strain = em;
            states(idx).global_stress   = Qbar * eg;
            states(idx).material_stress = Q * em;
        end
    end
end
