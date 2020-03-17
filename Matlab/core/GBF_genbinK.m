% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

function bf = GBF_genbinK(NodeInd, y, alpha)

% function [bf] = GBF_genbinK(U, Lambda, NodeInd, type, alpha)
%
% GBF_genbinK computes a Kernel Schur multiplier based on a binary classification
%
% In:
%    NodeInd   = K vector - The indices of the K basis nodes
%    y         = labels of the binary classification
%    alpha     = shape parameter of the binary kernel
%
% Out:
%    bf        = NxK matrix with Schur multipliers of binary kernel

N = length(y);
K = length(NodeInd);

% Initialize variables

bf = zeros(N,K);
idxpos = find(y>=0);
idxneg = find(y<0);

% Compute the matrix of basis vectors

for i=1:K
        if (y(NodeInd(i)) >=0)
           bf(idxpos,i)= ones(length(idxpos),1);
           bf(idxneg,i)= alpha*ones(length(idxneg),1);
        elseif (y(NodeInd(i)) < 0)
           bf(idxneg,i)= ones(length(idxneg),1);
           bf(idxpos,i)= alpha*ones(length(idxpos),1);  
        end
end
  
return