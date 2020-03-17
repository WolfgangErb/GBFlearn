% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

function [edges,A] = GBF_sim(nodes,alpha)

% Calculates edges and adjacency matrix for point clouds in order to
% generate a NN-graph G. 

% INPUT:  
% nodes        : Nodes of point cloud
% alpha        : Distance parameter for similarity graph
%
% OUTPUT:    
% edges        : The edges of the similarity graph
% A            : The adjacency matrix of the similarity graph

N = size(nodes,1);

% number of edges in a complete graph
M = N*(N-1)/2;

% determine the similarity weight of each edge
edges = zeros(M,2);
A = zeros(N,N);
iter = 0;

for i = 1 : N
    for j = i+1 : N
        iter = iter + 1;
        edges(iter,:) = [i j];
        A(i,j) = exp(-alpha*norm( nodes(i,:) - nodes(j,:), 2 ).^2);
        A(j,i) = exp(-alpha*norm( nodes(i,:) - nodes(j,:), 2 ).^2);
    end
end

return