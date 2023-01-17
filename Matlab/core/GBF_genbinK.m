% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

function bf = GBF_genbinK(idxW, y, alpha)

% function bf = GBF_genbinK(idxW, y, alpha)
%
% GBF_genbinK computes a Kernel Schur multiplier based on a binary classification
%
% Input:
%    idxW      = K vector - The indices of the K basis nodes
%    y         = labels of the binary classification
%    alpha     = shape parameter of the binary kernel
%
% Output:
%    bf        = NxK matrix with Schur multipliers of binary kernel

N = length(y);
K = length(idxW);

% Initialize variables

bf = zeros(N,K);
idxpos = find(y>=0);
idxneg = find(y<0);

% Compute the matrix of basis vectors

for i=1:K
        if (y(idxW(i)) >=0)
           bf(idxpos,i)= ones(length(idxpos),1);
           bf(idxneg,i)= alpha*ones(length(idxneg),1);
        elseif (y(idxW(i)) < 0)
           bf(idxneg,i)= ones(length(idxneg),1);
           bf(idxpos,i)= alpha*ones(length(idxpos),1);  
        end
end
  
return