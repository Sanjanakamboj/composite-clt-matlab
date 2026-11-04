function [eps0, kappa] = solve_deformation(ABD, N, M)
%SOLVE_DEFORMATION Mid-plane strain and curvature from applied resultants.
%   [eps0; kappa] = ABD \ [N; M]. Errors on a singular / ill-conditioned ABD.
    c = cond(ABD);
    if ~isfinite(c) || c > 1e14
        error('solve_deformation:illConditioned', ...
              'ABD matrix is singular or ill-conditioned (cond = %.3e).', c);
    end
    d = ABD \ [N(:); M(:)];
    eps0  = d(1:3);
    kappa = d(4:6);
end
