function [Z_final,A_final,B_final, F_final] = drjcc(X, k1, k2, W, options)
%%%%The model----------------
% min{Z,A,B,F}w1*||X'-Z'A||^2+alpha1||A||1+w2*||A-BF||^2+alpha2*Tr(FLF'), 
% s.t. B>=0,F>=0.

% You only need to provide the above five inputs.
% Notation:
% X ... (n x m) data matrix 
%       n  ... number of  genes (feature size)   
%       m  ... number of cell sample               
% k1 ... number of feature for dimension reduction
% k2 ... number of features for clustering in NMF(the number of cluster)
% W ... weight matrix of the affinity graph 
%
% options ... Structure holding all settings
%%%
% Coder Wenming Wu Email: wenmingwu55 at 163.com
%
%   version  --December/2019, revision,--March/2020

%   maxIter The maximum number of iterations
%   tol1,tol2 Iteration error
%   bSuccess
%     [n,m]=size(X);
    n = size(X,1);
    m = size(X,2);
    %%%%%%%%%%%%------Three initialization methods
    %%%------Random initialization
%     Z0= abs(rand(n,k1));
%     A0 = abs(rand(k1,m));
%     B0= abs(rand(k1,k2));
%     F0 = abs(rand(k2,m));

%%%------Initialization by SVD
    [U,vV,D] = svds(X,k1);
    Z0 = abs(U*sqrt(vV));
    A0 = abs(sqrt(vV)*D');
    
    [U,vV,D] = svds(A0,k2);
    B0 = abs(U*sqrt(vV));
    F0 = abs(sqrt(vV)*D');
    
%%%------Initialization by NMF
%     [Z0,A0] = snmf(X,k1,100);
%     [B0,F0] = snmf(A0,k2,100);
    
    
    D = diag(sum(W,2));
    L = D-W;
    
    Z = Z0;A = A0;
    B = B0;F = F0;
    E = A;
    T = zeros(k1,m);
    iter = 0; 
    converged = 0;    
    maxIter = 100;  
    tol1 = 1e-5;tol2 = 1e-5;
    tryNo = 0;
    w1 = 1;w2 = 1;
    
    %%%%%%%===========Update variables Z,A,B,F by iteration================

      while ~converged  && iter < maxIter   
%     tmp_T = cputime;
%     tryNo = tryNo+1;
%     nIter = 0;
%     maxErr = 1;
%     nStepTrial = 0;
    iter = iter + 1;
    derta = 5e+1;

    alpha1 = 0.5; alpha2 = 0.5;
        %%%========Update variables Z==========
    Zwk = Z;
    Z = Z.*((X*A')./(Z*A*A'+eps));
    Zkl = Z;
    %%%========Update variables A==========
    I = eye(k1);
    VV1 = w1*Z'*X+w2*B*F+derta*E+T;
    VV2 = (w1*Z'*Z+w2*I+derta*I)*A;
    Awk = A;
    A = A.*(VV1./VV2+eps);
    Akl = A;
    Ewk = E;
    E = soft(A-T/derta,alpha1/derta);
    Ekl = E;
    T = T+1.618*derta*(E-A);
    
    %%%%%Update variables B and F, where Bk and Fk are the variables at the k-th iteration, 
    %%%%% and Bkl and Fkl are the variables at the k+1-th iteration.==================
    %%%========Update variables B==========
    Bwk = B;
    B = B.*(A*F')./(B*F*F'+eps);
    B = max(B,0);
    Bkl = B;
    %%%========Update variables F==========
        
    Fwk = F;
    FF1 = w2*B'*A+alpha2*F*W;
    FF2 = w2*B'*B*F+alpha2*F*D;
    F = F.*(FF1)./(FF2);
    F = max(F,0);
    Fkl = F;

    temp = max ([norm(Zkl-Zwk,'fro'),norm(Akl-Awk,'fro'),norm(Bkl-Bwk,'fro'),norm(Fkl-Fwk,'fro')]);
    %     temp = muu*temp/norm(V,2);
    temp =temp/max([norm(X,'fro')]);
    
    if temp < tol1 
        converged = 1;
    end
    
%         disp([' µü´ú´ÎÊý ' num2str(iter) ' temp ' num2str(temp) ]);

    t1(iter)=temp;
    

end
    Z_final = Z; %%% Z_final  is finally Z
    A_final = A; %%% A_final  is finally A  
    B_final = B; %%% B_final  is finally B
    F_final = F; %%% F_final  is finally F
        
end
function[y] = soft( x, T )
  if sum( abs(T(:)) )==0
       y = x;
  else
       y = max( abs(x) - T, 0);
       y = sign(x).*y;
   end
end    