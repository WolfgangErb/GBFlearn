% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

function A = GBF_adjmat(edges,N)

% This function calculates the adjacency matrix of an unweighted 
% undirected graph from the information of the edges

% INPUT:    
% edges        : The edges of the graph G
% N            : Number of vertices of the graph G
%
% OUTPUT:  
% A            : Adjacency matrix of the graph G

A = zeros(N,N);
A(sub2ind([N N],edges(:,1),edges(:,2))) = 1;
A(sub2ind([N N],edges(:,2),edges(:,1))) = 1;

return