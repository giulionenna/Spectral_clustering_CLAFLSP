%% --- TASK 6-7-8---
M=3;
idx_circle={}; %cluster labels for circle dataset
idx_spiral={}; %same for spiral


for i = 1:length(K)
    idx_circle{i} = kmeans(U_circle{i}, M);
    idx_spiral{i} = kmeans(U_spiral{i}, M);
end
