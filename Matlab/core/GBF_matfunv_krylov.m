% GBFlearn: a toolbox for graph signal interpolation
% and classification with graph basis functions (GBFs)
% (C) W. Erb 15.01.2023

function yy = GBF_matfunv_krylov(A,vv,phi,m,tol,method)

% function yy = GBF_matfunv_krylov(A,vv,phi,m,tol,method)
%
% Computes the block vector yy = phi(A) vv for a matrix function phi 
% using an iterative block Lanczos procedure
%
% In:
%    A         = nxn symmetric matrix
%    vv        = nxk block vector consisting of orthogonal unit vectors
%    phi       = positive function on the interval [0,Lambda]
%    m         = number of iterations
%    tol       = tolerance
%    method    = applied block method 
%                (classical: 'cbl', global: 'gbl', sequential: 'sbl')
%
% Out:
%    yy        = nxk block vector

if ~exist('method','var')
     method = 'gbl';
end

[n,k] = size(vv);

switch method
    
    case 'cbl' %classical block Lanczos method
    
    V = zeros(n,k*(m+1));
    H = zeros(k*(m+1),k*(m+1));

    V(:,1:k) = vv;

    for j=1:m
        w = A*V(:,(j-1)*k+1:j*k);
        for p = 1:k
            for i=max(1,(j-2)*k+p):j*k+p-1
                H(i,(j-1)*k+p) = w(:,p)'*V(:,i);
                w(:,p)      = w(:,p) - H(i,(j-1)*k+p)*V(:,i);
            end
            H(j*k+p,(j-1)*k+p) = norm(w(:,p));
            V(:,j*k+p) = w(:,p)/H(j*k+p,(j-1)*k+p);
        end
    
        e1 = zeros(j*k,k); e1(1:k,1:k) = eye(k);
        ej = zeros(j*k,k); ej((j-1)*k+1:j*k,1:k) = eye(k);
    
        phiH      = matfun(H(1:j*k,1:j*k),phi);
        u         = phiH*e1;
        %MM = H(1:j*k,1:j*k);
        %u         = (MM+0.01*eye(j*k))\((MM+0.01*eye(j*k))\e1);
        resnorm   = norm(H(j*k+1:(j+1)*k,(j-1)*k+1:j*k)* (ej'*u));
    
        if (j >= 10) && (resnorm<=tol)
             break
        end

    end

    yy = V(:,1:j*k)*u;    
    
    case 'gbl'  %global block Lanczos method

    V = zeros(n*k,m+1);
    H = zeros(m+1,m+1);

    w = vv(:);
    beta = norm(w);
    V(:,1) = w/beta;

    for j=1:m
        psi = reshape(V(:,j),n,k);
        w = A*psi; w = w(:);
        H(j,j) = w'*V(:,j);
        for i=max(1,j-1):j
            w      = w - H(i,j)*V(:,i);
        end
        H(j+1,j) = norm(w); H(j,j+1) = H(j+1,j);
        V(:,j+1) = w/H(j+1,j);
    
        e1 = zeros(j,1); e1(1) = 1;
        ej = zeros(j,1); ej(j) = 1;
   
        u         = matfun(H(1:j,1:j),phi)*e1;
        resnorm   = abs(H(j+1,j)* (ej'*u));
    
        if (j >= 10) && (resnorm<=tol)
              break
        end

    end
    yy = V(:,1:j)*(beta*u);
    yy = reshape(yy,n,k);

    
    case 'sbl' %sequential block Lanczos method
    
    yy = zeros(n,k);
   
    for p = 1:k 
        
        V = zeros(n,m+1);
        beta = norm(vv(:,p));
        psi = vv(:,p)/beta;
        psi0 = zeros(n,1);
        V(:,1) = psi;
        H = zeros(m+1,m+1);
        
        for j=1:m       
            w = A*psi;
               
            H(j,j) = w'*psi;           
            w      = w - H(j,j)*psi;
            if (j > 1)
                w     = w - H(j-1,j)*psi0;
            end
            psi0 = psi;           
            H(j+1,j) = norm(w); H(j,j+1) = H(j+1,j);
            psi = w/H(j+1,j);
            V(:,j+1) = psi;
    
            e1 = zeros(j,1); e1(1) = 1;
            ej = zeros(j,1); ej(j) = 1;
   
            u         = matfun(H(1:j,1:j),phi)*e1;
            resnorm     = abs(H(j+1,j)* (ej'*u));
            
            if (j >= 10) && (resnorm<=tol)
                break
            end
        end
                  
        yy(:,p) = V(:,1:j)*u*beta;
        
    end

end

end

function B = matfun(A,phi)
    [U,Lambda] = GBF_spectrum(A,'ascend');
    f = phi(Lambda);
    B = U*diag(f)*U';
end

