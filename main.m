clear
close all
clc

eigs_true = true;
deflation = 'wiel';
print_fig = false; %Set TRUE for Printing figures to file
%% ---- LOADING DATA -----
load("data\Circle.mat", "X");
circle_X = X;
load("data\Spiral.mat", "X");
spiral_X = X;
clear X;

fig(1) = figure;
scatter(circle_X(:,1), circle_X(:,2), 10, 'filled');
set(fig(1),'PaperSize',[14 11]);
if print_fig==true
    print(fig(1), 'Latex\pictures\circle_scatterplot.pdf', '-dpdf')
end

fig(1) = figure;
scatter(spiral_X(:,1), spiral_X(:,2), 10, 'filled');
set(fig(1),'PaperSize',[14 11]);
if print_fig == true
    print(fig(1), 'Latex\pictures\spiral_scatterplot.pdf', '-dpdf')
end



%% --- TASK 1 ---
% --- BUILDING K-NN GRAPH ---
K=[10,20,40];   %Different value of K to test
G_circle={};    %Array containing Graph objects for circle dataset
W_circle={};    %Array containing adjacency matrices for circle dataset
G_spiral={};    
W_spiral={};

for i=1:length(K)
    [G_circle{i} , W_circle{i}]= knn_graph(circle_X, K(i), 1);
    [G_spiral{i}, W_spiral{i}] = knn_graph(spiral_X, K(i), 1);

    fig(2) = figure;
    plot(G_circle{i}, 'XData', circle_X(:,1), 'YData', circle_X(:,2))
    set(fig(2),'PaperSize',[14 11]);
    if print_fig == true
    print(fig(2), ['Latex\pictures\circle_KNN_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

    fig(2) = figure;
    plot(G_spiral{i}, 'XData', spiral_X(:,1), 'YData', spiral_X(:,2))
    set(fig(2),'PaperSize',[14 11]);
    if print_fig == true
        print(fig(2), ['Latex\pictures\spiral_KNN_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

end


%% --- TASK 2 ---

for i=1:length(K)
    [L_circle{i}, D_circle{i}] = graph_laplacian_norm(W_circle{i});
    [L_spiral{i}, D_spiral{i}] = graph_laplacian_norm(W_spiral{i});

    fig(3) = figure;
    spy(L_circle{i})
    set(fig(3),'PaperSize',[14 11]);
    if print_fig == true
        print(fig(3), ['Latex\pictures\circle_LaplacianSpy_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

    fig(3) = figure;
    spy(L_spiral{i})
    set(fig(3),'PaperSize',[14 11]);
    if print_fig == true
        print(fig(3), ['Latex\pictures\spiral_LaplacianSpy_K' int2str(K(i)) '.pdf'], '-dpdf')
    end
end
%% --- TASK 3-4-5 ---
n_eigs=6;
U_circle={};
U_spiral={};
eigs_circle={};
eigs_spiral={};

for i=1:length(K) %for each value of K tested
    if eigs_true
        [U_circle{i}, eigs_circle{i}] = eigs(L_circle{i}, n_eigs, 'smallestabs'); %compute the n_eigs smallest eigenvalues for circle
        [U_spiral{i}, eigs_spiral{i}] = eigs(L_spiral{i}, n_eigs, 'smallestabs'); %same for spiral
    else
        [eigs_circle{i}, U_circle{i}] = inverse_power_method_deflation(L_circle{i}, n_eigs, 1e-12, 1e4, deflation);
        [eigs_spiral{i}, U_spiral{i}] = inverse_power_method_deflation(L_spiral{i}, n_eigs, 1e-12, 1e4, deflation);
    end
end

for i = 1:length(K)
    %simple code to generate a matrix where in each column the eigenvalues
    %are listed for each value of K tested
    eigs_circle_tot(:,i) = diag(eigs_circle{i}); 
    eigs_spiral_tot(:,i) = diag(eigs_spiral{i});
end
M=3; %using the first 3 eigenvectors
for i = 1:length(K)
    %trim the matrix U to the first M columns
    U_circle{i} = U_circle{i}(:, 1:M);
    U_spiral{i} = U_spiral{i}(:, 1:M);
end




%% --- TASK 6-7-8---
M=3;
idx_circle={}; %cluster labels for circle dataset
idx_spiral={}; %same for spiral


for i = 1:length(K)
    idx_circle{i} = kmeans(U_circle{i}, M);
    idx_spiral{i} = kmeans(U_spiral{i}, M);
end

for i=1:length(K)

    fig(4) = figure;
    scatter(circle_X(:,1), circle_X(:,2), 10, idx_circle{i}, 'filled');
    set(fig(4),'PaperSize',[14 11]);
    if print_fig == true
    print(fig(4), ['Latex\pictures\circle_SpectralClustering_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

    fig(4) = figure;
    scatter(spiral_X(:,1), spiral_X(:,2), 10, idx_spiral{i}, 'filled');
    set(fig(4),'PaperSize',[14 11]);
    if print_fig == true
        print(fig(4), ['Latex\pictures\spiral_SpectralClustering_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

end

%% --- TASK 9 kmeans ---

%-- Kmeans no spectral---
idx_circle_kmeans = kmeans(circle_X, M);
idx_spiral_kmeans = kmeans(spiral_X(:,1:2), M);

fig(5) = figure;
scatter(circle_X(:,1), circle_X(:,2), 10, idx_circle_kmeans, 'filled');
set(fig(5),'PaperSize',[14 11]);
if print_fig == true
print(fig(5), ['Latex\pictures\circle_Kmeans_nospectral.pdf'], '-dpdf')
end
fig(5) = figure;
scatter(spiral_X(:,1), spiral_X(:,2), 10, idx_spiral_kmeans, 'filled');
set(fig(5),'PaperSize',[14 11]);
if print_fig == true
print(fig(5), ['Latex\pictures\spiral_Kmeans_nospectral.pdf'], '-dpdf')
end

%% --- Task 9 dbscan ---

idx_circle_dbscan = dbscan(circle_X, 0.6, 5);
idx_spiral_dbscan = dbscan(spiral_X(:,1:2), 0.3,5, 'Distance', 'mahalanobis');

fig(6) = figure;
scatter(circle_X(:,1), circle_X(:,2), 10, idx_circle_dbscan, 'filled');
set(fig(6),'PaperSize',[14 11]);
if print_fig == true
print(fig(6), ['Latex\pictures\circle_dbscan_nospectral.pdf'], '-dpdf')
end
fig(6) = figure;
scatter(spiral_X(:,1), spiral_X(:,2), 10, idx_spiral_dbscan, 'filled');
set(fig(6),'PaperSize',[14 11]);
if print_fig == true
print(fig(6), ['Latex\pictures\spiral_dbscan_nospectral.pdf'], '-dpdf')
end





%% test code

save('test_data', '-mat')

