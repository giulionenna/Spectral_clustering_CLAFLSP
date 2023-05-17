function [D, U] = power_method_deflation(A, K, tol, maxiter, factorization)
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
debug = 1;

%if issymmetric(A)
    [lambda_min, v_min] = inverse_power_method(A, tol, maxiter, 'lu');
%else
%    disp('Matrix is not symmetric')
%end

if(lambda_min<0)
    mu = abs(lambda_min)+tol;
else
    mu=0;
end
A = A+mu*eye(N); %shift spectrum so that A is spd

% Compute the new matrix on which the power method will be applied

for eig_iter = 1:K %for each eigenvalue to be computed
    [lambda, v] = power_method(A, tol, maxiter); %compute the greatest eigenvector of A
    D(eig_iter, eig_iter) = lambda-mu; %save eigenvalue of A to be shifted
    U(:, eig_iter) = v; %save eigenvector of A 

    % Deflate A
    A = A- lambda* v*v';
end

end

function [lambda, v] = inverse_power_method(A, tol, maxiter, factorization)  
    
    N = length(A);
    A = A+tol*eye(N);
    v_0 = rand(N,1);
    v_0 = v_0/norm(v_0);
    if strcmp(factorization, 'chol')
        [R, flag] = chol(A);
    else
       [L,U] = lu(A); 
    end
    v_old = v_0 / norm(v_0);
    lambda_old = 0;
    lambda_new = 2;
    k=0;
    while k<maxiter && abs((lambda_old - lambda_new)/lambda_old)>tol
        lambda_old = lambda_new;
        if strcmp(factorization, 'chol')
            z = R'\v_old;
            v_new = R\z;
        else
            z = L\v_old;
            v_new = U\z;
        end
        lambda_new = v_old'*v_new;
        v_old = v_new/norm(v_new);
        k = k+1;
    end
    lambda = 1/lambda_new - tol;
    v = v_new/norm(v_new);
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
