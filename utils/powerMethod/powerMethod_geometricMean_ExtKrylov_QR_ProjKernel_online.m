function [v, lambda, i] = powerMethod_geometricMean_ExtKrylov_QR_ProjKernel_online(A, B, M, numEigenVectors, LA, LB, powerMethod_numIter, krylov_numIter, krylov_TOL, powerMethod_TOL, MequalA, krylov_solverMaxIter, KernelMatrix)
% function [v, lambda, i] = powerMethod_geometricMean_ExtKrylov_QR_ProjKernel_online(A, B, M, numEigenVectors, LA, LB, numIter, numIterF, TOL_F, TOL_powerMethod, MequalA, solverMaxIter)
% [v, lambda] = powerMethod_geometricMean_ExtKrylov_QR_v1_pcg(A, B, M, LA, UA, numIter, numIterF, TOL_F, TOL_powerMethod)
% Power Method for the geometric mean of matrices A,B, i.e. A#B
% We consider the expression of A#B as A( A^-1 * B )^(1/2), leading to the
% following system of equations:
% For a given y we look for x such that
% A#Bx = y iff 
% A( A^-1 * B )^(1/2) x = y iff
% x = (A^-1 B)^(-1/2) A^-1 y
% IMPORTANT: in this version we aim to find simultaneously k eigenvectors,
% aided by QR decomposition

% This function returns the k eigenvectors correspoding to the smallest
% k eigenvalues of (A#B), based on the power method through a function of
% matrices and QR decomposition

KernelMatrixIsEmpty = false;

if size(KernelMatrix,2) > 1
    [KernelMatrixOrth, ~] = licols(KernelMatrix);
elseif size(KernelMatrix,2) == 1
    KernelMatrixOrth = KernelMatrix/norm(KernelMatrix);
elseif isempty(KernelMatrix)
    KernelMatrixOrth = [];
    KernelMatrixIsEmpty = true;
end

if nargin < 12
    krylov_solverMaxIter = [];
end

delta       = 1;
symm        = 0;
matrixPower = -1/2; % this is fixed, as for the moment we focuse only on the geometric mean

n           = size(A,1);
xkplus1     = zeros(n, numEigenVectors);

xk          = rand(n,numEigenVectors);
[xk,rk]     = qr(xk,0);
rk          = 1./diag(rk);


for i = 1:powerMethod_numIter
       
    for j = 1:numEigenVectors
        [y,~]         = pcg(A,xk(:,j), [], krylov_solverMaxIter, LA, LA');
        y             = real(y);
        xkplus1(:,j)  = gen_kpik_pcg_diagShift(A, B, M, y, matrixPower, delta, krylov_TOL, krylov_numIter, symm, LA, LB, MequalA, krylov_solverMaxIter);
    end
    
    xkplus1           = real(xkplus1);

    if ~KernelMatrixIsEmpty
        xkplus1 = xkplus1 - KernelMatrixOrth*(KernelMatrixOrth'*xkplus1);
    end
    
    [xkplus1,rkplus1] = qr(xkplus1,0);
    rkplus1           = 1./diag(rkplus1);
    
    stopCriterion     = sqrt(sum((1-abs(diag(xkplus1'*xk))).^2));
    xk                = xkplus1;
    rk                = rkplus1;

    if max(stopCriterion) < powerMethod_TOL
        break
    end

end
1;
v = xk;

lambda = rk;
1;
if i >= powerMethod_numIter
    warning('No convergence')
end
1;

