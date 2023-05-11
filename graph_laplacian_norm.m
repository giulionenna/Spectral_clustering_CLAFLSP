function [L, D] = graph_laplacian_norm(W)
    D = sparse(diag(sum(W).^(-1/2)));
    L = speye(size(W))-D*W*D;    
end