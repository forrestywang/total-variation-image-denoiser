% FormMatrix returns the FDM matrix A from vector u with weight alpha.
function A = FormMatrix(u, alpha)
    % Variables
    n = sqrt(length(u));
    h = 1 / (n + 1);
    c = - alpha / h^2;
    u = reshape(u, n, n);
    
    % Calculations
    % Compute neighbouring values
    % A_NE, A_SW are not needed for our stencil
    A_C = u;
    A = paddata(u, [n + 2, n + 2], Side = "both");
    A_N = A(1:end-2, 2:end-1);
    A_S = A(3:end, 2:end-1);
    A_E = A(2:end-1, 3:end);
    A_W = A(2:end-1, 1:end-2);
    A_NW = A(1:end-2, 1:end-2);
    A_SE = A(3:end, 3:end);
    
    % Compute stencil values
    AW = stencil(A_C - A_W, A_C - A_S, A_C - A_W, A_NW - A_W, h, c);
    AE = stencil(A_E - A_C, A_E - A_SE, A_E - A_C, A_N - A_C, h, c);
    AS = stencil(A_C - A_W, A_C - A_S, A_SE - A_S, A_C - A_S, h, c);
    AN = stencil(A_E - A_C, A_N - A_C, A_N - A_NW, A_N - A_C, h, c);
    AC = 1 - AW - AE - AS - AN;
    
    % Undo spdiags shifting
    AW = circshift(AW, -n);
    AN = circshift(AN, -1);
    AE = circshift(AE, n);
    AS = circshift(AS, 1);
    
    % Apply Dirichlet BCs
    AN(n:n:end) = 0;
    AS(n+1:n:end) = 0;
    
    % Output
    A = spdiags([AW AN AC AS AE], [-n -1 0 1 n], n^2, n^2);
end

% stencil is a helper function that helps compute the entries in the finite
% element stencil.
function A = stencil(A1, A2, A3, A4, h, c)
    % β > 0 is a constant parameter to circumvent numerical issues, such as
    % division by zero
    beta = 10^(-6);
    
    % Calculations
    mat1 = c ./ ( 2 * sqrt( (A1 / h).^2 + (A2 / h).^2 + beta ) );
    mat2 = c ./ ( 2 * sqrt( (A3 / h).^2 + (A4 / h).^2 + beta ) );
    
    % Output
    A = mat1 + mat2;
    A = A(:);
end