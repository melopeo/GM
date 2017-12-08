function [krylovOpts, eigsOpts] = process_inputs(krylovOpts, eigsOpts)

% Parameters for Krylov Method
if ~isfield(krylovOpts, 'solverMaxIter')
    krylovOpts.solverMaxIter = 200;
end

if ~isfield(krylovOpts, 'numIter')
    krylovOpts.numIter = 5;
end

if ~isfield(krylovOpts, 'TOL')
    krylovOpts.TOL = 1.e-5;
end

% Parameters for Power Method
if ~isfield(eigsOpts, 'numIter')
    eigsOpts.numIter = 1000;
end

if ~isfield(eigsOpts, 'TOL')
    eigsOpts.TOL = 1.e-6;
end
