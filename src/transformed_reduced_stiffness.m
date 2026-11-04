function Qbar = transformed_reduced_stiffness(mat, theta_deg)
%TRANSFORMED_REDUCED_STIFFNESS Qbar in global x-y axes for a ply at theta (deg).
%   Standard Jones/Reddy closed form, engineering shear strain convention.
    Q = reduced_stiffness(mat);
    Q11 = Q(1,1); Q12 = Q(1,2); Q22 = Q(2,2); Q66 = Q(3,3);
    th = theta_deg * pi / 180;
    m = cos(th); n = sin(th);
    m2 = m^2; n2 = n^2; m4 = m2^2; n4 = n2^2;
    Q11b = Q11*m4 + 2*(Q12 + 2*Q66)*m2*n2 + Q22*n4;
    Q22b = Q11*n4 + 2*(Q12 + 2*Q66)*m2*n2 + Q22*m4;
    Q12b = (Q11 + Q22 - 4*Q66)*m2*n2 + Q12*(m4 + n4);
    Q66b = (Q11 + Q22 - 2*Q12 - 2*Q66)*m2*n2 + Q66*(m4 + n4);
    Q16b = (Q11 - Q12 - 2*Q66)*m^3*n - (Q22 - Q12 - 2*Q66)*m*n^3;
    Q26b = (Q11 - Q12 - 2*Q66)*m*n^3 - (Q22 - Q12 - 2*Q66)*m^3*n;
    Qbar = [Q11b Q12b Q16b; Q12b Q22b Q26b; Q16b Q26b Q66b];
end
