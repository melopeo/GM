function y = gen_kpik_pcg_diagShift(X, Y, M, v, matrixPower, delta, tol, max_iter, symm, LX, LY, MequalX, solverMaxIter)
% function y = gen_kpik (X, Y, M, v, f, delta, tol, max_iter, symm)
%
% Compute an M-orthogonal basis of the generalized Krylov plus
% inverted Krylov subspace
% K_k(Q, P, v) = {v, Q \ (P v), P \ (Q v), ..., (Q \ P)^(k + 1) v}
% well as the generalized extended Arnoldi decomposition of P and Q.
% The value of the argument symm can be either 0 (use Arnoldi
% algorithm) or 1 (use Lanczos-like recursion).

% This is a modified version of the original shared by Massimiliano Fasi:
% - Weighted geometric mean of large-scale matrices: numerical analysis and algorithms 
% - http://amslaurea.unibo.it/8274/

% Computing the weighted geometric mean of two large-scale matrices and its inverse times a vector
% - http://eprints.ma.man.ac.uk/2474/

n = size(v, 1);
n_reorth = 2;

% Compute the maximum number of steps
m_max = min(floor (n / 2), max_iter);
m_max2 = 2 * m_max;
m_max22 = m_max2 + 2;

if (symm == 1)
    k_max = 4;
else
    k_max = m_max2;
end

sqrt_vMv = sqrt (v' * M * v); 

v1 = v;
[v2,~] = pcg(Y, X * v, [], solverMaxIter, LY, LY');
Vp = [v1 v2];

VpMVp    = Vp' * M * Vp;
tryAgain = true;
while tryAgain
    try
        Vp            = Vp / chol(VpMVp);
        tryAgain      = false;
    catch
        [~,eigMatrix] = eig(VpMVp,'vector');
        eigMatrix     = eigMatrix(:);
        minEig        = min(eigMatrix);
        minEig        = minEig*(minEig<0);
        diagShift     = 1.e-6 - minEig;
        VpMVp         = VpMVp + diagShift*eye(2);
    end
end
% try
%     Vp = Vp / chol((Vp' * M * Vp));
% catch
%     % sometimes the product Vp' * M * Vp is not positive definite. 
%     % Force it adding a small diagonal shift
%     [~,eigMatrix] = eig(Vp' * M * Vp,'vector');
%     eigMatrix     = eigMatrix(:);
%     minEig        = min(eigMatrix);
%     diagShift     = 1.e-6 - minEig;
%     Vp            = Vp / chol(Vp' * M * Vp + diagShift*eye(2)); 
% end
V(:,1 : 2) = Vp;

T = zeros(m_max22, m_max2);
D = zeros(n, m_max);

% choose the correct term for the computation of T
MmVm = zeros(n, m_max2);
if MequalX 
    MmVm(:, 1 : 2) = Y * V(:, 1 : 2);
else
    [vAux1, ~] = pcg(X, Y * V(:, 1), [], solverMaxIter, LX, LX');
    [vAux2, ~] = pcg(X, Y * V(:, 2), [], solverMaxIter, LX, LX');    
    MmVm(:, 1 : 2) = M*[vAux1 vAux2];
end
T(1 : 2, 1 : 2) = V(:, 1 : 2)' * MmVm(:, 1:2);

% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
for j = 1 : m_max
    j2   = j * 2;
    j22  = j2 + 2;
    j2m1 = j2 - 1;
    j21  = j2 + 1;

    [Vp(:, 1),~] = pcg(X, Y * V(:, j2m1), [], solverMaxIter, LX, LX'); 
    [Vp(:, 2),~] = pcg(Y, X * V(:, j2), [], solverMaxIter, LY, LY'); 

    % Extend orthogonal base
    k_min = max(1, j - k_max);
    for l = 1 : n_reorth
        for kk = k_min : j
            k1 = (kk - 1) * 2 + 1;
            k2 = kk * 2;
            coef = V(:, k1 : k2)' * M * Vp;
            Vp = Vp - V(:, k1 : k2) * coef;
        end
    end

    
    VpMVp    = Vp' * M * Vp;
    tryAgain = true;
    while tryAgain
        try
            Vp            = Vp / chol(VpMVp);
            tryAgain      = false;
        catch
            [~,eigMatrix] = eig(VpMVp,'vector');
            eigMatrix     = eigMatrix(:);
            minEig        = min(eigMatrix);
            minEig        = minEig*(minEig<0);
            diagShift     = 1.e-6 - minEig;
            VpMVp         = VpMVp + diagShift*eye(2);
        end
    end
%     try
%         Vp = Vp / chol((Vp' * M * Vp));
%     catch
%         % sometimes the product Vp' * M * Vp is not positive definite. 
%         % Force it adding a small diagonal shift
%         [~,eigMatrix] = eig(Vp' * M * Vp,'vector');
%         eigMatrix     = eigMatrix(:);
%         minEig        = min(eigMatrix);
%         diagShift     = 1.e-6 - minEig;
%         Vp            = Vp / chol(Vp' * M * Vp + diagShift*eye(2)); 
%     end 
    V(:, j21 : j22) = Vp;

    % Incremental update of T
    if (j ~= m_max)
        if MequalX %if (M == X)
            Vc = Y * V(:, j21:j22);
        else
            [vAux1, ~] = pcg(X, Y * V(:, j21), [], solverMaxIter, LX, LX');
            [vAux2, ~] = pcg(X, Y * V(:, j22), [], solverMaxIter, LX, LX');
            Vc = M * [vAux1 vAux2];
        end
        T(1 : j2, j21 : j22)    = V(:, 1 : j2)' * Vc;
        T(j21 : j22, 1 : j2)    = V(:, j21 : j22)' * MmVm(:, 1:j2);
        T(j21 : j22, j21 : j22) = V(:, j21 : j22)' * Vc;
        MmVm(:, j21 : j22)      = Vc;
    end

    [Vtemp,dtemp] = eig(T(1 : j2, 1 : j2),'vector');
    fT     = Vtemp*diag(dtemp.^matrixPower)*Vtemp';

    D(:, j) = V(:, 1 : j2) * (fT(:, 1) * sqrt_vMv); % Pedro

    % Estimate the residual during last step and return current u
    if (j >= delta + 1)
        dt = norm (D(:, j) - D(:, j - delta)) / norm (D(:, j - ...
        delta));
        if (abs (dt / (1 - dt)) <= tol)
            y = D(:, j);
            return
        end
    end
end

% If maximum number of iterations is reached without convergence, return
% the latest approximation
% warning('Convergence not reached')
y = D (:, end);