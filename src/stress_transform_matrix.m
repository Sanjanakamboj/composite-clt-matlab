function T = stress_transform_matrix(theta_deg)
%STRESS_TRANSFORM_MATRIX T(theta): rotates a GLOBAL stress vector into MATERIAL axes.
    th = theta_deg * pi / 180;
    m = cos(th); n = sin(th);
    T = [ m*m,  n*n,  2*m*n;
          n*n,  m*m, -2*m*n;
         -m*n,  m*n,  m*m - n*n];
end
