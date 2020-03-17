% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

function [s,c,Kf] = GBF_itpGBF(bf, NodeInd, y , delta)

% function [s,c,Kf] = GBF_itpGBF(bf, NodeInd, y, delta)
%
% Computes the GBF interpolant s, which
% fulfills the interpolation property
%            s(i) = y(i),
% for all sampling nodes i in NodeInd. 
%
% In:
%    bf        = Nxk matrix - the k basis vectors
%    NodeInd   = k vector - The indices of the k interpolation nodes
%    y         = k vector - The interpolation values at the k nodes
%    delta     = augmentation parameter for c.p.d. interpolation
%
% Out:
%    s         = N-vector - The GBF-interpolant
%    c         = k+M-vector - The coefficients (M=0 or M=1 describes 
%                wheter we solve with p.d. kernel or c.p.d kernel)
%    Kf        = interpolation matrix

if ~exist('delta','var')
    delta = 0;
end    

[N,k] = size(bf);

% Initialize variables

Kf = zeros(k,k);

% Generate matrix for the interpolation

for i=1:k
    Kf(i,:) = bf(NodeInd(i),:);
end

switch delta
    case 0     % standard p.d. kernel interpolation
        
    c  = Kf\y;
    s  = bf*c;
        
    case 'tps' % c.p.d kernel interpolation a la thin-plate spline
        
    Kf(k+1,1:k) = 1;
    Kf(1:k,k+1) = 1;
    y  = [y;0];
    
    c  = Kf\y;
    s  = [bf ones(N,1)]*c; 
    
    otherwise  % c.p.d kernel interpolation with augmented kernel
        
    Kf = Kf + ones(k)*delta;
    c  = Kf\y;
    s  = bf*c + delta*sum(c)*ones(N,1); 
             
end

return