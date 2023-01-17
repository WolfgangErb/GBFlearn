% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

function [phi,idxWsub,ysub] = GBF_genPUM(edges,cluster,domain,idxW,y)

% File:    GBF_genPUM.m
%
% Goal:    Computes PU based on given graph clusters and domains
%
% Use:     [phi,idxWsub,ysub] = GBF_genPUM(edges,cluster,domain,idxW,y)
%
% Input:
%           edges       = edges of the graph
%           cluster     = JJ clusters of the graph 
%           domain      = JJ augmented domains of the graph
%           idxW        = indices of the interpolation nodes
%           y           = interpolation values
%
% Output:
%           phi         = PU for subdomains
%           idxWsub     = indices of interpolation nodes in subdomains
%           ysub        = interpolation values in the subdomains

N = max(max(edges));
JJ = length(cluster);
phi = zeros(N,JJ);

for k = 1:JJ
    phi(cluster{k},k) = 1;
    idxWsub{k} = []; ysub{k} = [];  
    for i = 1:length(idxW)
        idx = find(domain{k} == idxW(i));
        if ~isempty(idx)
            idxWsub{k} = [idxWsub{k};idx];
            ysub{k} = [ysub{k};y(i)];
        end
    end
end

return
    