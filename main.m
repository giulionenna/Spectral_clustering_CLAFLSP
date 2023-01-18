clear
close all
clc

% ---- LOADING DATA -----
load("data\Circle.mat", "X");
circle_X = X;
load("data\Spiral.mat", "X");
spiral_X = X;
clear X;

fig(1) = figure
scatter(circle_X(:,1), circle_X(:,2), 10, 'filled');
set(fig(1),'PaperSize',[14 11]);
print(fig(1), 'Latex\pictures\circle_scatterplot.pdf', '-dpdf')

fig(2) = figure
scatter(spiral_X(:,1), spiral_X(:,2), 10, 'filled');
set(fig(2),'PaperSize',[14 11]);
print(fig(2), 'Latex\pictures\spiral_scatterplot.pdf', '-dpdf')




% --- TASK 1 ---
% --- BUILDING K-NN GRAPH ---
K = 10;
[G_circle , W_circle]= knn_graph(circle_X, K, 1);
[G_spiral, W_spiral] = knn_graph(spiral_X, K, 1);

figure(2)
plot(G_circle, 'XData', circle_X(:,1), 'YData', circle_X(:,2))
% figure(2)
%plot(G_spiral, 'XData', spiral_X(:,1), 'YData', spiral_X(:,2))

%--- TASK 2 ---
[L_circle, D_circle] = graph_laplacian(W_circle);
[L_spiral, D_spiral] = graph_laplacian(W_spiral);

%--- TASK 3-4-5 ---
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

%--- TASK 6-7---
idx_circle = kmeans(U_circle, M);
idx_spiral = kmeans(U_spiral, M);

%--- Task 8---
figure(3)
scatter(circle_X(:,1), circle_X(:,2), 10, idx_circle, 'filled')

figure(4)
scatter(spiral_X(:,1), spiral_X(:,2), 10, idx_spiral, 'filled')


