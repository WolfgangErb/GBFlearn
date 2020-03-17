% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

function [edges,A] = GBF_NN(nodes,r)

% Calculates edges and adjacency matrix for point clouds in order to
% generate a NN-graph G. 

% INPUT:  
% nodes        : Nodes of point cloud
% r            : Radius for generation of NN-edges
%
% OUTPUT:    
% edges        : The edges of the graph
% A            : The adjacency matrix of the graph

N = size(nodes,1);

% number of edges in a complete graph
M = N*(N-1)/2;

% determine the neighbors of each node
edges = zeros(M,2);
iter = 0;

for i = 1 : N
    for j = i+1 : N
        % if the distance between two nodes is less than r
        % then they are neighbors.
        if norm( nodes(i,:) - nodes(j,:), 2 ) <= r
            iter = iter + 1;
            edges(iter,:) = [i j];
        end
    end
end

% remove zeros in the indices
edges = edges(1:iter,:);

% construct the adjacency matrix for the local outputs
A = GBF_adjmat(edges,N);

return