

load("test_data.mat")
close all
L = L_circle{1};

[u_true, D_true] = eigs(L, 6, 'smallestabs');
d_true = diag(D_true);
[D_wiel, u_wiel] = inverse_power_method_deflation(L, 6, 1e-12, 1e4, 'wiel');
d_wiel = diag(D_wiel);
[D_naive, u_naive] = inverse_power_method_deflation(L, 6, 1e-12, 1e4, 'naive');
d_naive = diag(D_naive);

for i=1:6
    figure(i)
    plot([u_true(:,i), u_wiel(:,i), u_naive(:,i)])
end

