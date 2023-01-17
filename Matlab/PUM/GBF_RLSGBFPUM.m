% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

function s = GBF_RLSGBFPUM(LL,idxdomain,phi,idxWdomain,ydomain,type,alpha,lambda)
%
% File:    GBF_RLSGBFPUM.m
%
% Goal:    Computes a global approximant of the graph signal based on local regularized least squares 
%          GBF approximants and a partition of unity on the graph domain
%
% Use:     s = GBF_RLSGBFPUM(LL,idxdomain,phi,idxWdomain,ydomain,type,alpha,lambda)
%
% Input:
%          LL         = NxN matrix - graph Laplacian
%          idxdomain  = JJ subdomains of the PU
%          phi        = Partition of Unity
%          idxWdomain = indices of sampling nodes on subdomains
%          ydomain    = sampling values on subdomains
%          type       = type of kernel for local RLS
%          alpha      = kernel parameters
%          lambda     = regularization parameter for local RLS
%
% Output:
%          s         = N-vector - The RLS-GBFPUM solution

N = size(LL,1);
JJ = length(idxdomain);
s = zeros(N,1);

% Generate global GBFPUM interpolant
for k = 1:JJ
    % Calculate local interpolant
    L = LL(idxdomain{k},idxdomain{k});
    bf = GBF_genGBF(L,idxWdomain{k},type,alpha);
    sdomain = GBF_RLSGBF(bf, idxWdomain{k}, ydomain{k},lambda);
    
    % Add up to global approximant
    sdomainext = zeros(N,1); 
    sdomainext(idxdomain{k}) = sdomain;
    s = s + sdomainext.*phi(:,k);
end

return
    