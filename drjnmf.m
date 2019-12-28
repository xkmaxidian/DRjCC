function [DT_final,V_final,B_final, F_final] = drjnmf(X, k1, k2, W, options)
%
% where
%   X
% Notation:
% X ... (mFea x nSmp) data matrix 
%       mFea  ... number of  genes (feature size)   
%       nSmp  ... number of cell sample               
% k1 ... number of feature for dimension reduction
% k2 ... number of features for clustering in NMF
% W ... weight matrix of the affinity graph 
%
% options ... Structure holding all settings
%
% You only need to provide the above five inputs.
% min{Z,A,B,F}||X'-Z'A||^2+alpha1||A||1+||A-BF||^2+alpha2*Tr(FLF'), 
% s.t. B>=0,F>=0.
% Coder Wenming Wu Email: wenmingwu55 at 163.com
%
%   version  --December/2019

%   nIter_final  最终迭代次数
%   elapse_final 最终运行时间
%   bSuccess
%   objhistory_final 最终目标函数值

%     [m,nSmp]=size(X);
%     U = abs(rand(mFea,k));
%     V = abs(rand(nSmp,k)); 



differror = 1e-5;%微分误差？
if isfield(options,'error')
    differror = options.error;
end

maxIter = [];%最大迭代次数？
if isfield(options, 'maxIter')
    maxIter = options.maxIter;
end

nRepeat =1;%重复次数？
if isfield(options,'nRepeat')
    nRepeat = options.nRepeat;
end

minIterOrig = 30;%最小迭代初始？
if isfield(options,'minIter')
    minIterOrig = options.minIter;
end
minIter = minIterOrig-1;

meanFitRatio = 0.1;%最小匹配率？
if isfield(options,'meanFitRatio')
    meanFitRatio = options.meanFitRatio;
end

alpha = 1;
if isfield(options,'alpha')
    alpha = options.alpha;
end
%%%%%%===========参数设置===========
Norm = 2;
NormV = 1;


if min(min(X)) < 0
    error('Input should be nonnegative!');
end

[m,nSmp]=size(X);
DT0= abs(rand(m,k1));
V0 = abs(rand(k1,nSmp));
[mFea,nSmp]=size(V0);

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

W = alpha*W;%权矩阵
DCol = full(sum(W,2));%行元素之和构成列向量DCol
D = spdiags(DCol,0,speye(size(W,1)));%构成对角阵D
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
    B0 = abs(rand(mFea,k2));
    F0 = abs(rand(k2,nSmp));
else
    nRepeat = 1;
end

[DT0,V0] = NormalizeUV(DT0, V0', NormV, Norm);V0=V0';
[B0,F0] = NormalizeUV(B0, F0', NormV, Norm);F0=F0';

DTk=DT0;Vk=V0;
Bk=B0;Fk=F0;
Ak=Vk;
Tk= zeros(k1,nSmp);
iter = 0; 
converged = 0;    
maxIter=200;  
tol1=1e-5;tol2=1e-5;
tryNo=0;
w1=1;w2=1;
  while ~converged  && iter < maxIter   
    tmp_T = cputime;
    tryNo = tryNo+1;
    nIter = 0;
    maxErr = 1;
    nStepTrial = 0;
    iter = iter + 1;
        derta =5e+1;

        %%%%%=============数据投影X=D'V,更新D，V==================
%           lambda1=5e-3;lambda2=5e-3;

        lambda1=norm(X,1)/norm(DTk,'fro');lambda2=norm(Vk,1)/norm(Fk,'fro');%%%%%%%Regularization parameter
%          lambda1=norm(X,1)/norm(DTk,'fro');;lambda2=norm(Vk,2)/norm(Fk,'fro');
        %%%========更新D==========
        DTkl=DTk.*((X*Vk')./(DTk*Vk*Vk'));
        %%%========更新V==========
        I=eye(k1);
        VV1=w1*DTkl'*X+w2*Bk*Fk+derta*Ak+Tk;
        VV2=(w1*DTkl'*DTkl+w2*I+derta*I)*Vk;
        Vkl=Vk.*(VV1./VV2);
        Akl=soft(Vkl-Tk/derta,lambda1/derta);
        Tkl=Tk+1.618*derta*(Akl-Vkl);
        %%%%%=============数据投影V=B*F,更新B，F==================
        Bkl=Bk.*(Vkl*Fk')./(Bk*Fk*Fk');
        for i=1:size(Bkl,1)
            for j=1:size(Bkl,2)
                if Bkl(i,j)<0
                   Bkl(i,j)=0;
                else
                    Bkl(i,j)=Bkl(i,j);
                end
            end
        end
        FF1=w2*Bkl'*Vkl+lambda2*Fk*W;
        FF2=w2*Bkl'*Bkl*Fk+lambda2*Fk*D;
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
    
    DTwk = DTk;
    Vwk = Vk;
    Awk=Ak;
    Bwk = Bk;
    Fwk = Fk;
    
    DTk=DTkl;
    Vk=Vkl;
    Bk=Bkl;
    Fk=Fkl;
%%%%%%%%%%Error
%   Er1(iter,:)=abs(mean(X - DTk*Vk))./norm(DTk*Vk,'fro');
%   Er2(iter,:)=abs(mean(Vk - Bk*Fk))./norm( Bk*Fk,'fro');
%   er1(iter)=mean(Er1(iter,:));
%   er2(iter)=mean(Er2(iter,:));
%   er=er1+er2;
  
  
    temp = max ([norm(DTkl-DTwk,'fro'),norm(Vkl-Vwk,'fro'),norm(Bkl-Bwk,'fro'),norm(Fkl-Fwk,'fro')]);
    %     temp = muu*temp/norm(V,2);
    temp =temp/max([norm(X,'fro')]);
    %     temp = max([(sqrt(L)*norm(ZK-Zkm1,'fro')),norm(WK-Wkm1,'fro'),norm(EK-Ekm1,'fro')]);
    %     temp = muu*temp/norm(Y,'fro');
    %     
    %%%%%%%%%%%%%%%%%%
    temp1 = max(norm( (X - DTk*Vk),'fro'),norm( (Vk - Bk*Fk),'fro'))/max(norm( Bk*Fk,'fro'),norm( DTk*Vk,'fro'));
    if temp1 < tol1 && temp < tol2
    converged = 1;
    end
%     disp(['temp1 ',num2str(temp1)]);
%     disp([' 迭代次数 ' num2str(iter) ' temp1 ' num2str(temp1) ' temp ' num2str(temp)]);
    t1(iter)=temp1;
    t2(iter)=temp;
%     elapse = cputime - tmp_T;

end
        DT_final = DTkl;
        V_final = Vkl;   
        B_final = Bkl;
        F_final = Fkl;

        [DT_final,V_final] = NormalizeUV(DT_final, V_final', NormV, Norm);V_final=V_final';
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