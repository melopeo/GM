function [eigvecs, eigvals] = power_method_for_geometric_mean_kpik(L, Q, numEigenvectorsToCompute, kernelMatrix, krylovOpts, eigsOpts, parallelExecution)

% Krylov method parameters
krylov_solverMaxIter = krylovOpts.solverMaxIter;
krylov_numIter       = krylovOpts.numIter;
krylov_TOL           = krylovOpts.TOL;

% Power Method parameters
powerMethod_numIter  = eigsOpts.numIter;
powerMethod_TOL      = eigsOpts.TOL;

MequalA              = true;
LQ                   = get_ichol(Q);
LL                   = get_ichol(L);

if parallelExecution
    [eigvecs, eigvals] = ...
        powerMethod_geometricMean_ExtKrylov_QR_ProjKernel_online_par(Q, L, Q, numEigenvectorsToCompute, LQ, LL, powerMethod_numIter, krylov_numIter, krylov_TOL, powerMethod_TOL, MequalA, krylov_solverMaxIter, kernelMatrix);
else
    [eigvecs, eigvals] = ...
        powerMethod_geometricMean_ExtKrylov_QR_ProjKernel_online(Q, L, Q, numEigenvectorsToCompute, LQ, LL, powerMethod_numIter, krylov_numIter, krylov_TOL, powerMethod_TOL, MequalA, krylov_solverMaxIter, kernelMatrix);
end
