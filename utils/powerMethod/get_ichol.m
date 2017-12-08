function L = get_ichol(A)
% L = get_ichol(A)
% This function calculates the Sparse Incomplete Cholesky factorization. In
% case the matrix is nearly singular, it adds a small diagonal shift,
% until the factorization is feasible.

tryCholesky = true;
while tryCholesky
    try 
        L           = ichol(A);
        tryCholesky = false;
    catch
        tryCholesky = true;
        A           = A + (1.e-6)*eye(size(A));
        warning('Cholesky with A is problematic')
    end
end