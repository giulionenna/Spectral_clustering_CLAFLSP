%% --- TASK 3-4-5 ---
n_eigs=6;
U_circle={};
U_spiral={};
eigs_circle={};
eigs_spiral={};

for i=1:length(K) %for each value of K tested
    [U_circle{i}, eigs_circle{i}] = eigs(L_circle{i}, n_eigs, 'smallestabs'); %compute the n_eigs smallest eigenvalues for circle
    [U_spiral{i}, eigs_spiral{i}] = eigs(L_spiral{i}, n_eigs, 'smallestabs'); %same for spiral
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