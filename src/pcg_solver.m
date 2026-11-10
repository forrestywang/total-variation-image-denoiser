% PCG computes the approximate solution to Ax = b using the Preconditioned 
% Conjugate Gradient method. This will stop when it reaches the maximum 
% number of iterations or when ||r||₂ ≤ tol ||b||₂.
function [x, iter] = pcg_solver(A, b, x_initial, L, U, maxiter, tol)
    % Iteration 1
    x = x_initial;
    r = b - A * x;
    z = U \ (L \ r);
    
    p = z;
    alpha = dot(z, r) / dot(p, A * p);
    
    % Iterations 2 and onward
    for iter = 2:maxiter
        % Update
        denom = dot(z, r);
        x = x + alpha * p;
        r = r - alpha * A * p;
    
        % Check
        if norm(r) <= tol * norm(b)
            iter = iter - 1;
            return;
        end
    
        % A-inner product Gram-Schmidt
        z = U \ (L \ r);
        num = dot(z, r);
    
        beta = num / denom;
        p = z + beta * p;
        alpha = num / dot(p, A * p);
    end
    
    % Final iteration
    x = x + alpha * p;
end