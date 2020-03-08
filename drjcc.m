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




differror = 1e-5;
if isfield(options,'error')
    differror = options.error;
end

maxIter = [];
if isfield(options, 'maxIter')
    maxIter = options.maxIter;
end

nRepeat =1;
if isfield(options,'nRepeat')
    nRepeat = options.nRepeat;
end

minIterOrig = 30;
if isfield(options,'minIter')
    minIterOrig = options.minIter;
end
minIter = minIterOrig-1;

meanFitRatio = 0.1;
if isfield(options,'meanFitRatio')
    meanFitRatio = options.meanFitRatio;
end

alpha = 1;
if isfield(options,'alpha')
    alpha = options.alpha;
end
%%%%%%===========Parameter settings===========
Norm = 2;
NormV = 1;


if min(min(X)) < 0
    error('Input should be nonnegative!');
end

[n,m]=size(X);
Z0= abs(rand(n,k1));
A0 = abs(rand(k1,m));
[mFea,nSmp]=size(A0);

%%%%%%=========== Weight matrix setting===========

if isfield(options,'weight') && strcmpi(options.weight,'NCW')
    feaSum = full(sum(X,2));
    D_half = (X'*feaSum).^.5;
    for i = 1:nSmp
        X(:,i) = X(:,i)/D_half(i);
    end
end

if isfield(options,'alpha_nSmp') && options.alpha_nSmp
    alpha = alpha*nSmp;    
end

W = alpha*W;% Weight matrix
DCol = full(sum(W,2));% Sum of row elements constitutes column vector DCol
D = spdiags(DCol,0,speye(size(W,1)));% Compose Diagonal D
L = D - W;
if isfield(options,'NormW') && options.NormW
    D_mhalf = DCol.^-.5;

    tmpD_mhalf = repmat(D_mhalf,1,nSmp);
    L = (tmpD_mhalf.*L).*tmpD_mhalf';
    clear D_mhalf tmpD_mhalf;

    L = max(L, L');
end

bSuccess.bSuccess = 1;

%%%%%%%===========initialization================

    
selectInit = 1;
if ~exist('U','var')
    B0 = abs(rand(k1,k2));
    F0 = abs(rand(k2,m));
else
    nRepeat = 1;
end

[Z0,A0] = NormalizeUV(Z0, A0', NormV, Norm);A0=A0';
[B0,F0] = NormalizeUV(B0, F0', NormV, Norm);F0=F0';

Zk=Z0;Ak=A0;
Bk=B0;Fk=F0;
Ek=Ak;
Tk= zeros(k1,nSmp);
iter = 0; 
converged = 0;    
maxIter=100;  
tol1=1e-5;tol2=1e-5;
tryNo=0;
w1=1;w2=1;

%%%%%%%===========Update variables Z,A,B,F by iteration================

  while ~converged  && iter < maxIter   
%     tmp_T = cputime;
%     tryNo = tryNo+1;
%     nIter = 0;
%     maxErr = 1;
%     nStepTrial = 0;
    iter = iter + 1;
        derta =5e+1;

%%%%% Update variables Z and A, where Zk and Ak are the variables at the k-th iteration, 
%%%%% and Zkl and Akl are the variables at the k+1-th iteration.==================
        alpha1=norm(X,1)/norm(Zk,'fro');alpha2=norm(Ak,1)/norm(Fk,'fro');%%%%%%%Regularization parameter
%          lambda1=norm(X,1)/norm(DTk,'fro');;lambda2=norm(Vk,2)/norm(Fk,'fro');
        %%%========Update variables Z==========
        Zkl=Zk.*((X*Ak')./(Zk*Ak*Ak'));
        %%%========Update variables A==========
        I=eye(k1);
        VV1=w1*Zkl'*X+w2*Bk*Fk+derta*Ek+Tk;
        VV2=(w1*Zkl'*Zkl+w2*I+derta*I)*Ak;
        Akl=Ak.*(VV1./VV2);
        Ekl=soft(Akl-Tk/derta,alpha1/derta);
        Tkl=Tk+1.618*derta*(Ekl-Akl);
%%%%%Update variables B and F, where Bk and Fk are the variables at the k-th iteration, 
%%%%% and Bkl and Fkl are the variables at the k+1-th iteration.==================
        %%%========Update variables B==========
        Bkl=Bk.*(Akl*Fk')./(Bk*Fk*Fk');
        for i=1:size(Bkl,1)
            for j=1:size(Bkl,2)
                if Bkl(i,j)<0
                   Bkl(i,j)=0;
                else
                    Bkl(i,j)=Bkl(i,j);
                end
            end
        end
        %%%========Update variables F==========
        FF1=w2*Bkl'*Akl+alpha2*Fk*W;
        FF2=w2*Bkl'*Bkl*Fk+alpha2*Fk*D;
        Fkl=Fk.*(FF1)./(FF2);
        for i=1:size(Fkl,1)
            for j=1:size(Fkl,2)
                if Fkl(i,j)<0
                   Fkl(i,j)=0;
                else
                    Fkl(i,j)=Fkl(i,j);
                end
            end
        end
        
%         [DTkl,Vkl] = NormalizeUV(DTkl, Vkl', NormV, Norm);Vkl=Vkl';
        [Bkl,Fkl] = NormalizeUV(Bkl, Fkl', NormV, Norm);Fkl=Fkl';
    
    Zwk = Zk;
    Awk = Ak;
    Ewk=Ek;
    Bwk = Bk;
    Fwk = Fk;
    
    Zk=Zkl;
    Ak=Akl;
    Bk=Bkl;
    Fk=Fkl;
%%%%%%%%%% Error
%   Er1(iter,:)=abs(mean(X - DTk*Vk))./norm(DTk*Vk,'fro');
%   Er2(iter,:)=abs(mean(Vk - Bk*Fk))./norm( Bk*Fk,'fro');
%   er1(iter)=mean(Er1(iter,:));
%   er2(iter)=mean(Er2(iter,:));
%   er=er1+er2;
  
  
    temp = max ([norm(Zkl-Zwk,'fro'),norm(Akl-Awk,'fro'),norm(Bkl-Bwk,'fro'),norm(Fkl-Fwk,'fro')]);
    %     temp = muu*temp/norm(V,2);
    temp =temp/max([norm(X,'fro')]);
    %     temp = max([(sqrt(L)*norm(ZK-Zkm1,'fro')),norm(WK-Wkm1,'fro'),norm(EK-Ekm1,'fro')]);
    %     temp = muu*temp/norm(Y,'fro');
    %     
    %%%%%%%%%%%%%%%%%%
    temp1 = max(norm( (X - Zk*Ak),'fro'),norm( (Ak - Bk*Fk),'fro'))/max(norm( Bk*Fk,'fro'),norm( Zk*Ak,'fro'));
    if temp1 < tol1 && temp < tol2
    converged = 1;
    end
%     disp(['temp1 ',num2str(temp1)]);
%     disp([' µü´ú´ÎÊý ' num2str(iter) ' temp1 ' num2str(temp1) ' temp ' num2str(temp)]);
    t1(iter)=temp1;
    t2(iter)=temp;
%     elapse = cputime - tmp_T;

  end
  
        Z_final = Zkl; %%% Z_final  is finally Z
        A_final = Akl; %%% A_final  is finally A  
        B_final = Bkl; %%% B_final  is finally B
        F_final = Fkl; %%% F_final  is finally F

        [Z_final,A_final] = NormalizeUV(Z_final, A_final', NormV, Norm);A_final=A_final';
        [B_final,F_final] = NormalizeUV(B_final, F_final', NormV, Norm);F_final=F_final';
%         t=1:iter;
%         figure
%         plot(t,er,'r-'),xlabel('Iteration times');ylabel('Error');
end

function[y] = soft( x, T )
  if sum( abs(T(:)) )==0
       y = x;
  else
       y = max( abs(x) - T, 0);
       y = sign(x).*y;
   end
end    