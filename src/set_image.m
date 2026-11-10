% set_image returns two n×n grayscale images, one clear and one noisy 
% image. n_s_ratio controls the ratio of noise to signal. If the seed is 
% not chosen, one will be randomly generated.
function [image, noisy] = set_image(n, n_s_ratio, seed)
    % Seed
    if nargin < 3
        rng("shuffle");
    else
        rng(seed, "twister");
    end

    % Evaluate exact solution at point P = [x, y]
    h = 1 / (n + 1);
    xvec = h:h:1-h;
    yvec = h:h:1-h;
    [x, y] = meshgrid(xvec, yvec);
    r1 = sqrt( (x - 0.5).^2 + (y - 0.5).^2 );
    r2 = sqrt( (x - 0.5).^2/2 + (y - 0.5).^2 );
    image = -0.5 * (r1 < 0.2) + (r2 < 0.3);
    
    % Add noise
    stddev = norm(image, "fro") * n_s_ratio / n;
    noisy = image + stddev * randn(n, n);
end