function [G, W] = knn_graph_ND(X, K, sigma)
%An adjacency matrix of a K-nn graph has, for each column and row, only K
%nonzero elements that correspond to the nearest neighbors of each point.

N = length(X);


W = sparse(N, N);
for i = 1:N %for each point
    %Compute the similarity for every other point
    sim = zeros(N,1);
    for j = 1:N
        sim(j) = exp(-(norm(X(i,:)-X(j,:), 2)^2)/(2*sigma^2));
    end
    sim(i) = 0; %set similarity with itself to 0
    [sim, idx] = maxk(sim, K); %Compute the K most similar points and its indices
    W(i, idx) = sim; %set values of the adjacency matrix
    W(idx, i) = sim; 
end
G = graph(W);