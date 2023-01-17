% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

function bf = GBF_genGBFeff(L, idxW, type, alpha, method)

% function bf = GBF_genGBFeff(L, idxW, type, alpha, method)
%
% Computes K generalized translates of the GBF for the nodes in idxW
% Uses block Krylov methods for the calculation of the GBFs
%
% In:
%    L         = NxN matrix - the sparse graph Laplacian
%    idxW      = K vector - The indices of the K interpolation nodes
%    type      = type of graph basis function
%    alpha     = additional computational parameters for the iterations
%    method    = iterative method for the calculation of the GBFs
%
% Out:
%    bf        = NxK matrix with the K basis function vectors

N = size(L,1);
K = length(idxW);

%Initialize variables

bf = sparse(N,K);

% Compute the matrix of basis vectors

bf(sub2ind([N K],idxW,(1:K)')) = 1;

% Generate matrix for the interpolation

switch type
    
case 'varspline'      %Variation spline kernel

   alphastd = [2,0.001,100,1e-8];
   alpha = [alpha,alphastd(length(alpha)+1:length(alphastd))];
   
   if ~exist('method','var')
       if (mod(alpha(1),1) == 0)
           method = 'direct';
       else
           method = 'gbl';
       end
   end

   
   switch method
       
       case 'cheb'         %Chebyshev method based on interpolation
  
           fun = @(x) (alpha(2)+x).^(-alpha(1));
           bf = GBF_matfunv_cheb(L,bf,fun,alpha(3),alpha(4));
       
       case 'cheb2'        %squared Chebyshev method
  
           fun = @(x) (alpha(2)+x).^(-alpha(1)/2);
           bf = GBF_matfunv_cheb(L,bf,fun,alpha(3)/2,alpha(4));
           bf = GBF_matfunv_cheb(L,bf,fun,alpha(3)/2,alpha(4));
                  
       case 'cheblsq'      %Chebyshev method based on least squares 
      
          fun = @(x) (alpha(2)+x).^(-alpha(1));
          bf = GBF_matfunv_cheb(L,bf,fun,alpha(3),alpha(4),'lsq'); 
          
       case 'cbl'          %classical block Lanczos
           fun = @(x) (alpha(2)+x).^(-alpha(1));
           bf = GBF_matfunv_krylov(L,bf,fun,alpha(3),alpha(4),'cbl');   
          
       case 'gbl'          %global block Lanczos
           fun = @(x) (alpha(2)+x).^(-alpha(1));
           bf = GBF_matfunv_krylov(L,bf,fun,alpha(3),alpha(4),'gbl'); 
           
       case 'sbl'          %sequential block lanczos
           fun = @(x) (alpha(2)+x).^(-alpha(1));
           bf = GBF_matfunv_krylov(L,bf,fun,alpha(3),alpha(4),'sbl');
           
       case 'direct'       %direct inversion (for integer-valued s)
           if (sign(alpha(1))==1)  
              LS = (alpha(2)*speye(N)+L); 
              for i = 1:ceil(abs(alpha(1)))
                  bf = LS\bf;
              end
           elseif (sign(alpha(1))==-1)
              LS = (alpha(2)*speye(N)+L);  
              for i = 1:ceil(abs(alpha(1)))
                  bf = LS*bf;
              end
           end
   end
  
case 'diffusion'      %Diffusion kernel
   
   if ~exist('method','var')
       method = 'gbl';
   end
   alphastd = [80,50,1e-8];
   alpha = [alpha,alphastd(length(alpha)+1:length(alphastd))];
   
   switch method
       
       case 'cheb'         %Chebyshev method based on interpolation
      
          fun = @(x) exp(-alpha(1)*x);
          bf = GBF_matfunv_cheb(L,bf,fun,alpha(2),alpha(3));
       
       case 'cheb2'        %squared Chebyshev method
         
          fun = @(x) exp(-alpha(1)/2*x);
          bf = GBF_matfunv_cheb(L,bf,fun,alpha(2)/2,alpha(3));
          bf = GBF_matfunv_cheb(L,bf,fun,alpha(2)/2,alpha(3));
          
       case 'cheblsq'      %Chebyshev method based on least squares
      
          fun = @(x) exp(-alpha(1)*x);
          bf = GBF_matfunv_cheb(L,bf,fun,alpha(2),alpha(3),'lsq'); 
   
       case 'cbl'          %classical block Lanczos method
           
          fun = @(x) exp(-alpha(1)*x); 
          bf = GBF_matfunv_krylov(L,bf,fun,alpha(2),alpha(3),'cbl');
    
       case 'gbl'          %global block Lanczos method
           
          fun = @(x) exp(-alpha(1)*x); 
          bf = GBF_matfunv_krylov(L,bf,fun,alpha(2),alpha(3),'gbl');
          
       case 'sbl'          %sequential block Lanczos method
           
          fun = @(x) exp(-alpha(1)*x); 
          bf = GBF_matfunv_krylov(L,bf,fun,alpha(2),alpha(3),'sbl');
   end
   
case 'bandlimited'      %Bandlimited kernel
   
   if ~exist('method','var')
       method = 'gbl';
   end
   alphastd = [0.1,100,1e-8];
   alpha = [alpha,alphastd(length(alpha)+1:length(alphastd))];
   
   switch method
       
       case 'cheb'         %Chebyshev method based on interpolation
      
          fun = @(x)double(x<=alpha(1));
          bf = GBF_matfunv_cheb(L,bf,fun,alpha(2),alpha(3));
   
       case 'cheb2'        %squared Chebyshev method
         
          fun = @(x)double(x<=alpha(1));
          bf = GBF_matfunv_cheb(L,bf,fun,alpha(2)/2,alpha(3));
          bf = GBF_matfunv_cheb(L,bf,fun,alpha(2)/2,alpha(3));
          
       case 'cheblsq'      %Chebyshev method based on least squares
      
          fun = @(x)double(x<=alpha(1));
          bf = GBF_matfunv_cheb(L,bf,fun,alpha(2),alpha(3),'lsq'); 
           
       case 'cbl'          %classical block Lanczos method
           
          fun = @(x)double(x<=alpha(1));
          bf = GBF_matfunv_krylov(L,bf,fun,alpha(2),alpha(3),'cbl');
          
       case 'gbl'          %global block Lanczos method
           
          fun = @(x)double(x<=alpha(1)); 
          bf = GBF_matfunv_krylov(L,bf,fun,alpha(2),alpha(3),'gbl');
          
       case 'sbl'          %sequential block Lanczos method
           
          fun = @(x)double(x<=alpha(1));
          bf = GBF_matfunv_krylov(L,bf,fun,alpha(2),alpha(3),'sbl');
   end

  
case 'trivial'        %Trivial kernel
      
  % do nothing as bf is already calculated
      
end

  
return