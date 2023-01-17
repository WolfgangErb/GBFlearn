% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

function B = GBF_lanczosbidiag(L,n_lanczos,v)

% function B = GBF_lanczosbidiag(L,n_lanczos,v)
%
% Computes low-dimensional bidiagonal matrix B to approximate norm of L
%
% In:
%    L         = nxn matrix 
%    n_lanczos = dimension of bidiagonal matrix
%    v         = initial vector
%
% Out:
%    B         = bidiagonal matrix of dimension n_lanczos+1 x n_lanczos 

B = zeros(n_lanczos+1,n_lanczos);
d = zeros(size(v));

beta = norm(v); u = v/beta;
for i=1:n_lanczos
  rl = L'*u - beta*d;
  alph = norm(rl); d = rl/alph;
  B(i,i) = alph;
  p = L*d - alph*u;
  beta = norm(p); u = p/beta;
  B(i+1,i) = beta;
end

return