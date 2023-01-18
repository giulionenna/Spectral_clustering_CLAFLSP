function [W] = adj_matrix(X, sigma, tol)
    
    X = X(:, 1:2); % N x 2 matrix containing point coordinates
    N = length(X);

    b = nchoosek(N, 2); %number of all possible combinations between points
    k = 1;
    v = zeros(b, 3); %matrix containing all nonzero elements of the adjacency matrix
                     % and its relative positions

    for i = 1:N-1
        for j = i+1:N
            s = exp(- (norm(X(i,:)-X(j,:))^2)/(2*sigma^2));
            if s>tol
                v(k, 1) = s;
            end
            v(k, 2) = i;
            v(k, 3) = j;
            k = k+1;
        end
    end
    v(k, :) = [0, N, N];
    W = sparse(v(:,2), v(:,3), v(:,1));
end