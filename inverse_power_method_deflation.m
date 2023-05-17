function [D, U] = inverse_power_method_deflation(A, K, tol, maxiter, factorization )
%INVERSE_POWER_METHOD computes eigenvalues and eigenvectors of a given matrix using the inverse power method algorithm
%   A = SYIMMETRIC POSITIVE DEFINITE matrix
%   K = number of eigenvectors and eigenvalues to be computed
%   factorization = type of factorization to be used for solving the linear
%   system
%
N = length(A);
D = spalloc(K, K, K);
U = zeros([N, K]);
eig_iter = 1;

%compute smallest eigenvalue and relative eigenvector
[lambda, v] = inverse_power_method(A, tol, maxiter, factorization);
D(eig_iter, eig_iter) = lambda;
U(:,eig_iter) = v;

%compute the householder matrix P such that Pv=e1
u = U(:,eig_iter);
u(1) = sign(u(1))*norm(u);
P = eye(N) - (2/norm(u)^2)*(u*u');

%compute similar matrix to A (same eigenvalues)
B = P*A*P';
A_2 = B(2:end,2:end);
[lambda, v] = inverse_power_method(A_2, tol, maxiter, factorization);


eig_iter = eig_iter+1;
end

function [lambda, v] = inverse_power_method(A, tol, maxiter, factorization)
    N = length(A);
    v_0 = rand(N,1);
    if strcmp(factorization, 'chol')
        [R, flag] = chol(A);
    end
    v_old = v_0 / norm(v_0);
    lambda_old = 0;
    lambda_new = 2;
    k=0;
    while k<maxiter && abs((lambda_old - lambda_new)/lambda_old)>tol
        lambda_old = lambda_new;
        z = R'\v_old;
        v_new = R\z;
        lambda_new = v_old'*v_new;
        v_old = v_new/norm(v_new);
        k = k+1;
    end
    lambda = 1/lambda_new;
    v = v_new/norm(v_new);
end