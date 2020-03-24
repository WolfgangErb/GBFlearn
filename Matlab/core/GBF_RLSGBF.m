% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

function [s,c,Kf] = GBF_RLSGBF(bf, idxW, y , lambda)

% function [s,c,Kf] = GBF_RLSGBF(bf, idxW, y, lambda)
%
% GBF_RLSGBF computes the GBF regularized least squares solution s,
% minimizing the functional
%            1/K \sum_i (s(i) - y(i)) + \lambda \|s\|_(K),
% based on K sampling nodes in idxW. 
%
% In:
%    bf        = NxK matrix - the K graph basis vectors
%    idxW      = K vector - The indices of the K sampling nodes
%    y         = K vector - The sampling values at the K nodes
%    lambda    = regularization parameter for RLS
%
% Out:
%    s         = N-vector - The GBF-RLS solution
%    c         = K - The coefficients in the basis expansion
%    Kf        = interpolation matrix of the kernel

[~,K] = size(bf);

% Initialize variables

Kf = zeros(K,K);

% Generate matrix for the interpolation

for i=1:K
    Kf(i,:) = bf(idxW(i),:);
end

Kfreg = Kf + eye(K)*lambda*K;
c  = Kfreg\y;
s  = bf*c; 

return