function [L, D] = graph_laplacian(W)
    D = sparse(diag(sum(W)));
    L = D-W;
end