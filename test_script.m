rng(420)

A = rand(100);


A = A'*A;
d_true = eigs(A, 5, 'smallestabs')
d_max_true = eigs(A, 1)
d = inverse_power_method_deflation(A, 5, 1e-12, 100000000)