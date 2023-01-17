% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

function [idxsub] = GBF_domainaugment(edges,idxsub,RR)

% File:     GBF_domainaugment.m
%
% Goal:     Augments the graph cluster domains
%
% Use:      [idxsub] = GBF_domainaugment(edges,idxsub,RR)
%
% Input:
%           edges  = edges of the graph
%           idxsub = JJ subdomains of the graph 
%           RR     = extension radius for augmentation
%
% Output:
%           idxsub = augmentend subdomains of the graph

JJ = length(idxsub);
for r = 1:RR
   for j = 1:JJ
      idxsub1 = edges(ismember(edges(:,1),idxsub{j}),2);
      idxsub2 = edges(ismember(edges(:,2),idxsub{j}),1);
      idxsub{j} = unique([idxsub{j},idxsub1',idxsub2']);
   end
end
return