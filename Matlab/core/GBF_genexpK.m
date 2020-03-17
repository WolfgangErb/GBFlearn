% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

function bf = GBF_genexpK(NodeInd, y, alpha)

% function [bf] = GBF_genexpK(U, Lambda, NodeInd, type, alpha)
%
% GBF_genexpK computes a Kernel Schur multiplier based on an exponentially weighted
% similarity graph
%
% In:
%    NodeInd   = K vector - The indices of the K basis nodes
%    y         = type of graph basis function
%    alpha     = shape parameter of the exponential kernel
%
% Out:
%    bf        = NxK matrix with Schur multipliers of exponential kernel

N = size(y,1);
K = length(NodeInd);

% Initialize variables

bf = zeros(N,K);

% Compute the matrix of basis vectors

for i=1:K
      for j = 1:N
           bf(j,i)= exp(-alpha*norm(y(NodeInd(i),:)-y(j,:))^2);
      end
end
  
return