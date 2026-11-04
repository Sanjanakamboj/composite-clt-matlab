function [ABD, z, areal_mass] = laminate_abd(mat, angles_deg)
%LAMINATE_ABD Assemble the 6x6 ABD matrix of a laminate (plies listed bottom-to-top).
%   mat        : ply material struct (E1,E2,G12,nu12,density,thickness)
%   angles_deg : vector of ply angles, bottom to top
%   ABD        : [A B; B D], units N/m, N, N*m
%   z          : ply interface coordinates relative to the mid-plane (m)
%   areal_mass : kg/m^2
    nply = numel(angles_deg);
    t = mat.thickness;
    h = nply * t;
    z = -h/2 + (0:nply) * t;
    A = zeros(3); B = zeros(3); D = zeros(3);
    for k = 1:nply
        Qb = transformed_reduced_stiffness(mat, angles_deg(k));
        A = A + Qb * (z(k+1)   - z(k));
        B = B + Qb * (z(k+1)^2 - z(k)^2);
        D = D + Qb * (z(k+1)^3 - z(k)^3);
    end
    B = B / 2;  D = D / 3;
    ABD = [A B; B D];
    areal_mass = nply * mat.density * t;
end
