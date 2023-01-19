clear
close all
clc

print_fig = true; %Set TRUE for Printing figures to file
%% ---- LOADING DATA -----
load("data\Circle.mat", "X");
circle_X = X;
load("data\Spiral.mat", "X");
spiral_X = X;
clear X;

fig(1) = figure;
scatter(circle_X(:,1), circle_X(:,2), 10, 'filled');
set(fig(1),'PaperSize',[14 11]);
if print_fig
    print(fig(1), 'Latex\pictures\circle_scatterplot.pdf', '-dpdf')
end

fig(1) = figure;
scatter(spiral_X(:,1), spiral_X(:,2), 10, 'filled');
set(fig(1),'PaperSize',[14 11]);
if print_fig
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
    if print_fig
    print(fig(2), ['Latex\pictures\circle_KNN_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

    fig(2) = figure;
    plot(G_spiral{i}, 'XData', spiral_X(:,1), 'YData', spiral_X(:,2))
    set(fig(2),'PaperSize',[14 11]);
    if print_fig
    print(fig(2), ['Latex\pictures\spiral_KNN_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

end


%% --- TASK 2 ---
for i=1:length(K)
    [L_circle{i}, D_circle{i}] = graph_laplacian(W_circle{i});
    [L_spiral{i}, D_spiral{i}] = graph_laplacian(W_spiral{i});

    fig(3) = figure;
    spy(L_circle{i})
    set(fig(3),'PaperSize',[14 11]);
    if print_fig
    print(fig(3), ['Latex\pictures\circle_LaplacianSpy_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

    fig(3) = figure;
    spy(L_spiral{i})
    set(fig(3),'PaperSize',[14 11]);
    if print_fig
    print(fig(3), ['Latex\pictures\spiral_LaplacianSpy_K' int2str(K(i)) '.pdf'], '-dpdf')
    end
end
%% --- TASK 3-4-5 ---
M=3;
[U_circle, eigs_circle] = eigs(L_circle, 6, 'smallestabs');
%circle has clearly 2 connected components since the multiplicity of
%eigenvalue 0 is 2. Eigenvalues that are "almost 0" are relative to "almost
%connected components" Hence the number of cluster should be equal to the
%number of eigenvalues that are "almost 0"
[U_spiral, eigs_spiral] = eigs(L_spiral, 6, 'smallestabs');
%Same for spiral but here we only have one connected component and 3
%clusters.

U_circle = U_circle(:,1:M);
U_spiral = U_spiral(:,1:M);

%% --- TASK 6-7---
idx_circle = kmeans(U_circle, M);
idx_spiral = kmeans(U_spiral, M);

%% --- Task 8---
figure(3)
scatter(circle_X(:,1), circle_X(:,2), 10, idx_circle, 'filled')

figure(4)
scatter(spiral_X(:,1), spiral_X(:,2), 10, idx_spiral, 'filled')


