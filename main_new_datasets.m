clear
close all
clc

print_fig = true; %Set TRUE for Printing figures to file
%% ---- LOADING DATA -----
load("data\tetra.mat");
tetra_X = tetra{2:end,{'xAxis', 'yAxis', 'zAxis'}};
%tetra_X = tetra_X(2:end, :);
load("data\chainlink.mat");
chainlink_X = chainlink{2:end,{'xAxis', 'yAxis', 'zAxis'}};
%chainlink_X = chainlink_X(2:end, :);
load("data\atom.mat")
atom_X = atom{2:end,{'xAxis', 'yAxis', 'zAxis'}};
%atom_X = atom_X(2:end, :);

fig(1) = figure;
scatter3(tetra_X(:,1), tetra_X(:,2), tetra_X(:,3), 10, tetra{2:end,"Cluster_id"}, 'filled');
set(fig(1),'PaperSize',[14 11]);
if print_fig==true
    print(fig(1), 'Latex\pictures\tetra_scatterplot.pdf', '-dpdf')
end


fig(1) = figure;
scatter3(chainlink_X(:,1), chainlink_X(:,2), chainlink_X(:,3), 10, chainlink{2:end,"Cluster_id"}, 'filled');
set(fig(1),'PaperSize',[14 11]);
if print_fig==true
    print(fig(1), 'Latex\pictures\chainlink_scatterplot.pdf', '-dpdf')
end


fig(1) = figure;
scatter3(atom_X(:,1), atom_X(:,2), atom_X(:,3), 10, atom{2:end,"Cluster_id"}, 'filled');
set(fig(1),'PaperSize',[14 11]);
if print_fig==true
    print(fig(1), 'Latex\pictures\atom_scatterplot.pdf', '-dpdf')
end


%% --- TASK 1 ---
% --- BUILDING K-NN GRAPH ---
K=[5,10,20];       %Different value of K to test
G_chainlink={};    %Array containing Graph objects for chainlink dataset
W_chainlink={};    %Array containing adjacency matrices for chainlink dataset
G_tetra={};    
W_tetra={};
G_atom={};    
W_atom={};

for i=1:length(K)
    [G_tetra{i} , W_tetra{i}]= knn_graph_ND(tetra_X, K(i), 1);
    [G_chainlink{i}, W_chainlink{i}] = knn_graph_ND(chainlink_X, K(i), 1);
    [G_atom{i}, W_atom{i}] = knn_graph_ND(atom_X, K(i), 1);
    
    fig(2) = figure;
    plot(G_tetra{i}, 'XData', tetra_X(:,1), 'YData', tetra_X(:,2), 'ZData', tetra_X(:,3))
    set(fig(2),'PaperSize',[14 11]);
    if print_fig == true
    print(fig(2), ['Latex\pictures\tetra_KNN_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

    fig(2) = figure;
    plot(G_chainlink{i}, 'XData', chainlink_X(:,1), 'YData', chainlink_X(:,2), 'ZData', chainlink_X(:,3))
    set(fig(2),'PaperSize',[14 11]);
    if print_fig == true
    print(fig(2), ['Latex\pictures\chainlink_KNN_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

    fig(2) = figure;
    plot(G_atom{i}, 'XData', atom_X(:,1), 'YData', atom_X(:,2), 'ZData', atom_X(:,3))
    set(fig(2),'PaperSize',[14 11]);
    if print_fig == true
    print(fig(2), ['Latex\pictures\atom_KNN_K' int2str(K(i)) '.pdf'], '-dpdf')
    end
end


%% --- TASK 2 ---

for i=1:length(K)
    [L_tetra{i}, D_tetra{i}] = graph_laplacian(W_tetra{i});
    [L_chainlink{i}, D_chainlink{i}] = graph_laplacian(W_chainlink{i});
    [L_atom{i}, D_atom{i}] = graph_laplacian(W_atom{i});

    fig(3) = figure;
    spy(L_tetra{i})
    set(fig(3),'PaperSize',[14 11]);
    if print_fig == true
        print(fig(3), ['Latex\pictures\tetra_LaplacianSpy_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

    fig(3) = figure;
    spy(L_chainlink{i})
    set(fig(3),'PaperSize',[14 11]);
    if print_fig == true
        print(fig(3), ['Latex\pictures\chainlink_LaplacianSpy_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

    fig(3) = figure;
    spy(L_atom{i})
    set(fig(3),'PaperSize',[14 11]);
    if print_fig == true
        print(fig(3), ['Latex\pictures\chainlink_atom_K' int2str(K(i)) '.pdf'], '-dpdf')
    end
end
%% --- TASK 3-4-5 ---
n_eigs=10;
U_tetra={};
U_chainlink={};
U_atom={};
eigs_tetra={};
eigs_chainlink={};
eigs_atom={};
for i=1:length(K) %for each value of K tested
    [U_tetra{i}, eigs_tetra{i}] = eigs(L_tetra{i}, n_eigs, 'smallestabs'); %compute the n_eigs smallest eigenvalues for tetra
    [U_chainlink{i}, eigs_chainlink{i}] = eigs(L_chainlink{i}, n_eigs, 'smallestabs'); %same for chainlink
    [U_atom{i}, eigs_atom{i}] = eigs(L_atom{i}, n_eigs, 'smallestabs'); %same for atom
end

for i = 1:length(K)
    %simple code to generate a matrix where in each column the eigenvalues
    %are listed for each value of K tested
    eigs_tetra_tot(:,i) = diag(eigs_tetra{i}); 
    eigs_chainlink_tot(:,i) = diag(eigs_chainlink{i});
    eigs_atom_tot(:,i) = diag(eigs_atom{i});
end

for i = 1:length(K)
    %trim the matrix U to the first M columns
    U_tetra{i} = U_tetra{i}(:, 1:4); %using 4 for tetra
    U_chainlink{i} = U_chainlink{i}(:, 1:2); %using 2 for chainlink
    U_atom{i} = U_atom{i}(:, 1:2); %matrix is badly conditioned, computation of eigs is inaccurate
end




%% --- TASK 6-7-8---
M=3;
idx_tetra={}; %cluster labels for tetra dataset
idx_chainlink={}; %same for chainlink
idx_atom={};


for i = 1:length(K)
    idx_tetra{i} = kmeans(U_tetra{i}, 4);
    idx_chainlink{i} = kmeans(U_chainlink{i}, 2);
    idx_atom{i} = kmeans(U_atom{i}, 2);
end

for i=1:length(K)

    fig(4) = figure;
    scatter3(tetra_X(:,1), tetra_X(:,2), tetra_X(:,3), 10, idx_tetra{i}, 'filled');
    set(fig(4),'PaperSize',[14 11]);
    if print_fig == true
    print(fig(4), ['Latex\pictures\tetra_SpectralClustering_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

    fig(4) = figure;
    scatter3(chainlink_X(:,1), chainlink_X(:,2), chainlink_X(:,3), 10, idx_chainlink{i}, 'filled');
    set(fig(4),'PaperSize',[14 11]);
    if print_fig == true
        print(fig(4), ['Latex\pictures\chainlink_SpectralClustering_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

    fig(4) = figure;
    scatter3(atom_X(:,1), atom_X(:,2), atom_X(:,3), 10, idx_atom{i}, 'filled');
    set(fig(4),'PaperSize',[14 11]);
    if print_fig == true
    print(fig(4), ['Latex\pictures\atom_SpectralClustering_K' int2str(K(i)) '.pdf'], '-dpdf')
    end

end



