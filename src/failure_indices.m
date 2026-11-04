function [fi, mode] = failure_indices(sigma, S)
%FAILURE_INDICES Ply failure indices from material-axis stress [s1; s2; t12] (Pa).
%   S: strength struct with positive magnitudes Xt, Xc, Yt, Yc, S12.
%   fi.max_stress, fi.tsai_hill, fi.tsai_wu; mode = governing Max-Stress mode.
    s1 = sigma(1); s2 = sigma(2); t12 = sigma(3);

    % Maximum Stress
    if s1 >= 0, r1 = s1 / S.Xt; else, r1 = -s1 / S.Xc; end
    if s2 >= 0, r2 = s2 / S.Yt; else, r2 = -s2 / S.Yc; end
    r12 = abs(t12) / S.S12;
    [fi.max_stress, which] = max([r1 r2 r12]);
    labels = {'fiber_tension','fiber_compression','transverse_tension', ...
              'transverse_compression','shear'};
    if which == 1, mode = labels{1 + (s1 < 0)};
    elseif which == 2, mode = labels{3 + (s2 < 0)};
    else, mode = labels{5};
    end

    % Tsai-Hill (Azzi-Tsai-Hill, allowable picked by stress sign per axis)
    if s1 >= 0, X = S.Xt; else, X = S.Xc; end
    if s2 >= 0, Y = S.Yt; else, Y = S.Yc; end
    fi.tsai_hill = (s1/X)^2 - (s1*s2)/X^2 + (s2/Y)^2 + (t12/S.S12)^2;

    % Tsai-Wu, F12 = -0.5*sqrt(F11*F22) (conventional, not experimentally validated)
    F1  = 1/S.Xt - 1/S.Xc;   F2  = 1/S.Yt - 1/S.Yc;
    F11 = 1/(S.Xt*S.Xc);     F22 = 1/(S.Yt*S.Yc);   F66 = 1/S.S12^2;
    F12 = -0.5 * sqrt(F11 * F22);
    fi.tsai_wu = F1*s1 + F2*s2 + F11*s1^2 + F22*s2^2 + F66*t12^2 + 2*F12*s1*s2;
end
