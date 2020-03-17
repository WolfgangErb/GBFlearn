% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 01.03.2020

function [bf,pdf] = GBF_genGBF(U, Lambda, NodeInd, type, alpha, epsilon)

% function [bf,pdf] = GBF_genGBF(U, Lambda, NodeInd, type, alpha, epsilon)
%
% Computes K generalized translates of the GBF for the nodes in NodeInd
%
% In:
%    U         = NxN matrix - the Fourier transform matrix on G
%    Lambda    = Nx1 vector - the eigenvalues of the Laplacian
%    NodeInd   = K vector - The indices of the K interpolation nodes
%    type      = type of graph basis function
%    alpha     = additional shape parameter
%    epsilon   = additional shape parameter
%
% Out:
%    bf        = NxK matrix with the K basis function vectors
%    pdf       = the generating p.d. function f

N = length(Lambda);
K = length(NodeInd);

%Initialize variables

bf = zeros(N,K);
base = zeros(N,1);

% Generate matrix for the interpolation

switch type
    
  case 'varspline'      %Variation spline kernel
  if ~exist('epsilon','var')
      epsilon = 0;
  end    
    
  f = (epsilon+Lambda).^(alpha);
  f(abs(f)>=1e12) = 0;
  
  case 'diffusion'      %Diffusion kernel
  
  f = exp(alpha*Lambda);
  
  case 'polydecay'      %Kernel with polynomial decay
      
  f = (1:N).^(alpha); f = f';
  
  case 'bandlimited'    %Bandlimited kernel
     
  f = [ones(K,1);zeros(N-K,1)];
  
  case 'trivial'        %Trivial kernel
  
  f = ones(N,1);

end

 A = U*diag(f)*U';
 pdf = U*f;

% Compute the matrix of basis vectors

for i=1:K
    base(NodeInd(i))=1;
    bf(:,i)=A*base;
    base(NodeInd(i))=0;
end
  
return