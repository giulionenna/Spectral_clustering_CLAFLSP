%SAME AS TEST SCRIPT BUT USING THE NORMALIZED LAPLACIAN 
clear
close all
clc
load("test_data_norm.mat")
close all
print_fig = true
L = L_circle{1};

tic
[u_true, D_true] = eigs(L, 15, 'smallestabs');
time_eigs = toc
d_true = diag(D_true);
tic
[D_wiel, u_wiel, times_wiel] = inverse_power_method_deflation(L, 15, 1e-12, 1e4, 'wiel');
time_wiel = toc
d_wiel = diag(D_wiel);
tic
[D_naive, u_naive, times_naive] = inverse_power_method_deflation(L, 15, 1e-12, 1e4, 'naive');
time_naive = toc
d_naive = diag(D_naive);

for i=1:6
    fig(i) = figure;
    plot([u_true(:,i), u_wiel(:,i), u_naive(:,i)], 'LineWidth',1.2)
    legend('eigs', 'wiel deflation', 'naive deflation')
    set(fig(i), 'PaperSize', [14, 14]);
    if print_fig == true
        print(fig(i), ['Latex\pictures\ipmd_test\eigenvector_norm_', int2str(i), '.pdf'], '-dpdf')
    end
end

fig(7) = figure;
plot([d_true], 'LineStyle', 'none', 'Marker','o', 'MarkerSize', 10, 'LineWidth',2)
hold on
plot([d_naive], 'LineStyle', 'none','Marker', 'x',  'MarkerSize', 10, 'LineWidth',2)
hold on
plot(d_wiel, 'LineStyle', 'none','Marker','*', 'MarkerSize', 10, 'LineWidth',2)
legend({'eigs', 'naive deflation', 'wiel deflation'}, 'Location', 'northwest')
grid on
if print_fig == true
    set(fig(7), 'PaperSize', [14, 14]);
    print(fig(7), ['Latex\pictures\ipmd_test\eigenvalues_comp_norm.pdf'], '-dpdf')
end

fig(8) = figure;
plot([times_naive, times_wiel], 'LineWidth',2)
legend('naive deflation', 'wiel deflation')
if print_fig == true
    set(fig(8), 'PaperSize', [14, 14]);
        print(fig(8), ['Latex\pictures\ipmd_test\times_norm'], '-dpdf')
    end