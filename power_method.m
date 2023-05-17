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


