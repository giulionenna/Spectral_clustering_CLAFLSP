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