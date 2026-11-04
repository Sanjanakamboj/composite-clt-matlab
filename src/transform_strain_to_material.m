function eps_mat = transform_strain_to_material(eps_global, theta_deg)
%TRANSFORM_STRAIN_TO_MATERIAL Global engineering strain -> material engineering strain.
%   Uses the Reuter matrix R = diag(1,1,2): eps_12 = R*T*inv(R)*eps_xy.
    R    = diag([1 1 2]);
    Rinv = diag([1 1 0.5]);
    T = stress_transform_matrix(theta_deg);
    eps_mat = R * T * Rinv * eps_global(:);
end
