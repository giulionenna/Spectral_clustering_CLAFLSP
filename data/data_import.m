%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: C:\Users\giuli\Dropbox (Politecnico Di Torino Studenti)\POLITO 2022_2023\Computational Linear Algebra FLSP\Spectral_clustering_CLAFLSP\data\tetra.csv
%
% Auto-generated by MATLAB on 11-May-2023 12:19:27

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["id", "xAxis", "yAxis", "zAxis", "Cluster_id"];
opts.VariableTypes = ["double", "double", "double", "double", "categorical"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Cluster_id", "EmptyFieldRule", "auto");

% Import the data
tetra = readtable("C:\Users\giuli\Dropbox (Politecnico Di Torino Studenti)\POLITO 2022_2023\Computational Linear Algebra FLSP\Spectral_clustering_CLAFLSP\data\tetra.csv", opts);


%% Clear temporary variables
clear opts

save('tetra.mat', '-mat')