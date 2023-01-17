% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

function [bf,pdf] = GBF_genGBF(L, idxW, type, alpha)

% function [bf,pdf] = GBF_genGBF(L, idxW, type, alpha)
%
% Computes K generalized translates of the GBF for the nodes in idxW
% Uses the spectral decomposition of the graph Laplacian L
%
% In:
%    L         = NxN matrix - the graph Laplacian
%    idxW      = K vector - The indices of the K interpolation nodes
%    type      = type of graph basis function
%    alpha     = additional shape parameter
%
% Out:
%    bf        = NxK matrix with the K basis function vectors
%    pdf       = the generating p.d. function f

N = size(L,1);
K = length(idxW);

%Calculate Spectrum of Laplacian
[U,Lambda] = GBF_spectrum(L,'ascend');

%Initialize variables

bf = zeros(N,K);
base = zeros(N,1);

% Generate matrix for the interpolation

switch type
    
  case 'varspline'      %Variational spline kernel
  if length(alpha)==1
      alpha(2) = 0;
  end    
    
  f = (alpha(2)+Lambda).^(-alpha(1));
  f(abs(f)>=1e12) = 0;
  
  case 'diffusion'      %Diffusion kernel
  
  f = exp(-alpha(1)*Lambda);
  
  case 'polydecay'      %Kernel with polynomial decay
  if length(alpha)==1
      alpha(2) = 1;
  end      
      
  f = transpose(1+alpha(2)*(0:N-1)).^(-alpha(1));
  
  case 'bandlimited'    %Bandlimited kernel
     
  f = [ones(K,1);zeros(N-K,1)];
  
  case 'trivial'        %Trivial kernel
  
  f = ones(N,1);

end

 A = U*diag(f)*U';
 pdf = U*f;

% Compute the matrix of basis vectors

for i=1:K
    base(idxW(i))=1;
    bf(:,i)=A*base;
    base(idxW(i))=0;
end
  
return