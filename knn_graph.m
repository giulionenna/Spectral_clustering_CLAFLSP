function [G, W] = knn_graph(X, K, sigma)

% W full is the "full" adjacency matrix with nonzero elements only in the
% upper diagonal

%An adjacency matrix of a K-nn graph has, for each column and row, only K
%nonzero elements that correspond to the nearest neighbors of each point.

N = length(X);
X = X(:,1:2);

W = sparse(N, N);
for i = 1:N %for each point
    %Compute the similarity for every other point
    sim = exp(-((X(i, 1)- X(:, 1)).^2 + (X(i, 2)- X(:, 2)).^2)/(2*sigma^2));
    sim(i) = 0; %set similarity with itself to 0
    [sim, idx] = maxk(sim, K); %Compute the K most similar points and its indices
    W(i, idx) = sim; %set values of the adjacency matrix
    W(idx, i) = sim; 
    
end
G = graph(W);