% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

function [idxcluster,idxQ] = GBF_Jcenters_greedy(A,J,idxW,idxQin)

% File:    GBF_Jcenters_greedy.m
%
% Goal:    Performs a reduced J-center clustering based on greedy search: 
%          Finds clusters based on a reduced greedy procedure 
%          on the subset W (given by idxW) of the nodes
%
% Use:     [idxcluster,idxQ] = GBF_Jcenters_greedy(A,J,idxW,idxQin)
%
% Input:
%           A          = adjacency matrix
%           J          = number of clusters
%           idxW       = indices of subset on which the greedy search is performed
%           idxQin     = predefined input centers (optional)
%
% Output:
%           idxcluster = Jx1 cell vector with cluster assignments
%           idxQ       = Jx1 vector of cluster center indices

N = size(A,1);
D = zeros(J,N);

if (nargin == 3)
    idxQin = idxW(randi([1,length(idxW)],1,1));
end

sizeidxQin = length(idxQin);
idxQ = zeros(J,1); 
idxQ(1:sizeidxQin) = idxQin;

for j = 1:J-1
    D(j,:) = distcenter(A,N,idxQ(j));
    if j >= sizeidxQin 
       [~,indmax] = max(min(D(1:j,idxW),[],1));
       idxQ(j+1) = idxW(indmax(1));
    end
end
D(J,:) = distcenter(A,N,idxQ(J));

[~,cluster] = min(D,[],1);
for k = 1:J
    idxcluster{k} = find(cluster==k);
end

return

% distcenter calculates the graph geodesic distance 
% of the graph nodes to the centers indicated by idxQ
% (calculation based on the adjacency matrix A)

function dist = distcenter(A,N,idxQ)

  s = zeros(N,1); s(idxQ) = 1;
  dist = zeros(N,1);
  step = 1;
  while sum(s)~=N
      R = A*s;
      new = (R ~=0).*(s==0);
      dist = dist + step*new;
      s = s + new;
      step = step + 1;
  end 
return