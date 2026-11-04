function Q = reduced_stiffness(mat)
%REDUCED_STIFFNESS Plane-stress reduced stiffness Q in the ply 1-2 axes (Pa).
%   mat: struct with fields E1, E2, G12, nu12 (SI units).
%   Engineering shear strain convention: [s1;s2;t12] = Q*[e1;e2;g12].
    nu21  = mat.nu12 * mat.E2 / mat.E1;
    denom = 1 - mat.nu12 * nu21;
    Q11 = mat.E1 / denom;
    Q22 = mat.E2 / denom;
    Q12 = mat.nu12 * mat.E2 / denom;
    Q66 = mat.G12;
    Q = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66];
end
