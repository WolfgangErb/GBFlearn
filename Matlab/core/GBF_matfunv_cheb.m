% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

function s = GBF_matfunv_cheb(L,v,fun,m,Tol,method)

% function s = GBF_matfunv_cheb(L,v,fun,m,Tol,method)
%
% Computes the block vector s = fun(A) v for a matrix function fun 
% using a Chebyshev method
%
% In:
%    L         = nxn symmetric matrix
%    v         = nxk block vector consisting of orthogonal unit vectors
%    fun       = positive function on the interval [0,Lambda]
%    m         = number of iterations
%    Tol       = tolerance
%    method    = applied block method 
%                (interpolation: 'itp', least squares: 'lsq', non-negative: 'lsqnn')
%
% Out:
%    yy        = nxk block vector

if ~exist('method','var')
   method = 'itp'; 
end

N = size(L,1);

% Compute estimate of norm of L with Lanczos bidiagonalization
B = GBF_lanczosbidiag(L,10,v(:,1));
Lnorm = 1.1*norm(B);

range = [0,Lnorm];

% Calculate Chebyshev expansion coefficients for exponential
if strcmp(method,'itp') == 1
   coeff = chebcoeff(fun,m,range);
elseif strcmp(method,'lsq') == 1
   coeff = chebcoeffLSQ(fun,m,range,'lsq');
elseif strcmp(method,'lsqnn') == 1
   coeff = chebcoeffLSQ(fun,m,range,'lsqnn');
end
  
% Calculate Chebyshev approximation of exponential via recursion
L = (2*L - (range(2)+range(1))*speye(N))/(range(2) - range(1));

v_old = v;
v = L*v;
s = coeff(1)*v_old + coeff(2)* v;

for i = 2:m
    v_new = 2*L*v-v_old;
    s = s + coeff(i+1)*v_new;
    if norm(coeff(i+1)*v_new(:,1)) < Tol*norm(s(:,1)) 
        break; 
    end
    v_old = v;
    v = v_new;
end

% fprintf('Number of iteration steps: %2d \n',i);

end

function coeff = chebcoeff(fun,m,range)
  % Calculate Chebyshev interpolation coefficients of function fun
  zx = cos(transpose(linspace(0,1,m + 1))*pi);
  zx = (zx+1)*(range(2)-range(1))/2 + range(1);
  w = ones(m+1,1); w(1) = w(1)/2; w(m+1) = w(m+1)/2;

  f = fun(zx);
  gh = real(fft(f.*w,2*m));   
  gh = gh(1:m+1);
  coeff = gh.*w*2/m;
end

function coeff = chebcoeffLSQ(fun,m,range,par)
  % Calculate Chebyshev least squares coefficients of function fun
  
  if ~exist('par','var')
     par = 'lsq'; 
  end
  
  K = 2;         % multplication factor for number of nodes

  zx = cos(linspace(0,1,K*m + 1)'*pi);
  IM = cos(acos(zx)*(0:m));
  
  zx = (zx+1)*(range(2)-range(1))/2 + range(1);
  f = fun(zx);
  
  coeff = IM\f;
  
  if strcmp(par,'lsqnn') == 1
      K2 = 4;
      zx2 = cos(linspace(0,1,K2*m + 1)'*pi);
      IM2 = cos(acos(zx2)*(0:m));
      if (nnz(IM2*coeff < -eps) > 0)
         coeff = lsqlin(IM,f,-IM2,zeros(K2*m+1,1));
      end
  end

end