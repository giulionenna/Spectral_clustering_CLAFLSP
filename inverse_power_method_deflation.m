function [D, U] = inverse_power_method_deflation(A, K, tol, maxiter, deflation_type)
%INVERSE_POWER_METHOD computes eigenvalues and eigenvectors of a given matrix using the inverse power method algorithm
%   A = SYIMMETRIC POSITIVE DEFINITE matrix
%   K = number of eigenvectors and eigenvalues to be computed
%   factorization = type of factorization to be used for solving the linear
%   system
%
N = length(A);
D = spalloc(K, K, K);
U = zeros([N, K]);

% Compute the greatest eigenvalue of the matrix A
[lambda_max, v_max] = power_method(A, tol, maxiter);

% Compute the new matrix on which the power method will be applied
% (spectrum shift)
mu = lambda_max+1e3*tol;
B = mu*eye(N)-A;

if strcmp(deflation_type, 'naive')
    for eig_iter = 1:K %for each eigenvalue to be computed
        [lambda, v] = power_method(B, tol, maxiter); %compute the greatest eigenvector of B
        D(eig_iter, eig_iter) = mu-lambda; %save eigenvalue of A from the back of the spectrum
        U(:, eig_iter) = v; %save eigenvector of A (is the same of B)
        % Deflate B
        B = B - lambda* (v)*(v');
    end
elseif strcmp(deflation_type, 'wiel')
    % performing eigenvalues computation using Wielandt deflation
    for eig_iter = 1:K %for each eigenvalue to be computed
        [lambda, v] = power_method(B, tol, maxiter); %compute the greatest eigenvector of B
        D(eig_iter, eig_iter) = mu-lambda; %save eigenvalue of A from the back of the spectrum
        U(:, eig_iter) = v; %save eigenvector of A (is the same of B)
        % Deflate B
        [~, idx] = max(v);
        idx = min(idx);
        x = (1/lambda*v(idx))*(B(idx,:)');
        B = B - lambda* (v)*(x');
    end
    
else
    disp('Specify a valid deflation type')
end


end


function [lambda, v] = power_method(A, tol, maxiter)
    N = length(A);
    v_0 = rand(N,1);
    v_old = v_0 / norm(v_0);
    lambda_old = 0;
    lambda_new = 2;
    k=0;
    while k<maxiter && abs((lambda_old - lambda_new)/lambda_old)>tol
        lambda_old = lambda_new;
        v_new = A*v_old;
        lambda_new = v_old'*v_new;
        v_old = v_new/norm(v_new);
        k = k+1;
    end
    lambda = lambda_new;
    v = v_new/norm(v_new);
end