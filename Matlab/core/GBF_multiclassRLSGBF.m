% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

function sclass = GBF_multiclassRLSGBF(bf, idxW, y , lambda)

% function sclass = GBF_multiclassRLSGBF(bf, idxW, y, lambda)
%
% GBF_RLSGBF computes the GBF regularized least squares solution s,
% minimizing the functional
%            1/K \sum_i (s(i) - y(i)) + \lambda \|s\|_(K),
% based on K sampling nodes in idxW. 
% 
% GBF_multiclassRLSGBF uses a one-vs-rest strategy to reduce
% multiclass problems to binary classification problems with GBF_RLSGBF 
% as a basic solver.
% 
% In:
%    bf        = NxK matrix - the K graph basis vectors
%    idxW      = K vector - The indices of the K sampling nodes
%    y         = K vector - The sampling values at the K nodes
%    lambda    = regularization parameter for RLS
%
% Out:
%    sclass    = N-vector - The GBF-RLS multiclass solution

N = size(bf,1);

% Determine number of classes and initialize
classlabel = unique(y);
RLSvalue   = sum(abs(classlabel));
classN     = length(classlabel);            
s          = zeros(N,classN);
sclass     = zeros(N,1);

% Perform one-vs-rest strategy 
for i=1:classN
    classx = y;
    classx(classx==classlabel(i)) = RLSvalue;
    classx(classx~=RLSvalue) = -1;
    classx(classx==RLSvalue) = 1;
    s(:,i) = GBF_RLSGBF(bf,idxW,classx,lambda);
end

% Determine classes according to highest confidence score 
[~,index] = max(s,[],2);

for i = 1:classN
    sclass(index==i) = classlabel(i);
end


return