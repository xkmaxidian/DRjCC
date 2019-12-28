function [AMI_]=AMI(true_mem,mem)
%Program for calculating the Adjusted Mutual Information (AMI) between
%two clusterings.
%
% This code is a modified version of Nguyen Xuan Vinh's one available online.
%
% The modification includes the computation of AMI for contingency tables
% with empty rows or columns.
%
%--------------------------------------------------------------------------
%**Input: a contingency table T
%   OR
%        cluster label of the two clusterings in two vectors
%        eg: true_mem=[1 2 4 1 3 5]
%                 mem=[2 1 3 1 4 5]
%        Cluster labels are coded using positive integer. 
%**Output: AMI: adjusted mutual information  (AMI normalized by Sqrt(HA,HB))
%
%**Note: In a prevous published version, if you observed strange AMI results, eg. AMI>>1, 
%then it's likely that in these cases the expected MI was incorrectly calculated (the EMI is the sum
%of many tiny elements, each falling out the precision range of the computer).
%However, you'll likely see that in those cases, the upper bound for the EMI will be very
%tiny, and hence the AMI -> NMI (see [3]). It is recommended setting AMI=NMI in
%these cases, which is implemented in this version.
%--------------------------------------------------------------------------
%References: 
% [1] 'A Novel Approach for Automatic Number of Clusters Detection based on Consensus Clustering', 
%       N.X. Vinh, and Epps, J., in Procs. IEEE Int. Conf. on 
%       Bioinformatics and Bioengineering (Taipei, Taiwan), 2009.
% [2] 'Information Theoretic Measures for Clusterings Comparison: Is a
%	    Correction for Chance Necessary?', N.X. Vinh, Epps, J. and Bailey, J.,
%	    in Procs. the 26th International Conference on Machine Learning (ICML'09)
% [3] 'Information Theoretic Measures for Clusterings Comparison: Variants, Properties, 
%       Normalization and Correction for Chance', N.X. Vinh, Epps, J. and
%       Bailey, J., Journal of Machine Learning Research, 11(Oct), pages
%       2837-2854, 2010


if nargin==1
    T=true_mem; %contingency table pre-supplied
elseif nargin==2
    %build the contingency table from membership arrays
    r=max(true_mem);
    c=max(mem);

    %identify & removing the missing labels
    list_t=ismember(1:r,true_mem);
    list_m=ismember(1:c,mem);
    T=Contingency(true_mem,mem);
    T=T(list_t,list_m);
end

[r c]=size(T);
if (c == 1 || r == 1)
 error('Clusterings should have at least 2 clusters')
 return
end

N = sum(sum(T)); % total number of records

% update the true dimensions
a=sum(T,2)';
b=sum(T);

%calculating the Entropies
Ha=-(a(a ~= 0)/N)*log2(a(a ~= 0)/N)'; 
Hb=-(b(b ~= 0)/N)*log2(b(b ~= 0)/N)';

%calculate the MI (unadjusted)
MI=0;
for i=1:r
    for j=1:c
        if T(i,j)>0 MI=MI+T(i,j)*log2(T(i,j)*N/(a(i)*b(j)));end;
    end
end
MI=MI/N;

%-------------correcting for agreement by chance---------------------------
AB=a'*b;
bound=zeros(r,c);

E3=(AB/N^2).*log2(AB/N^2);
E3(isnan(E3)) = 0; % substitute 0log0=NaN with 0s

EPLNP=zeros(r,c);
log2Nij=log2([1:min(max(a),max(b))]/N);
for i=1:r
    for j=1:c
        nij=max(1,a(i)+b(j)-N);
        X=sort([nij N-a(i)-b(j)+nij]);
        if N-b(j)>X(2)
            nom=[[a(i)-nij+1:a(i)] [b(j)-nij+1:b(j)] [X(2)+1:N-b(j)]];
            dem=[[N-a(i)+1:N] [1:X(1)]];
        else
            nom=[[a(i)-nij+1:a(i)] [b(j)-nij+1:b(j)]];       
            dem=[[N-a(i)+1:N] [N-b(j)+1:X(2)] [1:X(1)]];
        end
        p1=prod(nom./dem)/N;                
        for nij=max(1,a(i)+b(j)-N):1:min(a(i), b(j))
            EPLNP(i,j)=EPLNP(i,j)+nij*log2Nij(nij)*p1;            
            p1=p1*(a(i)-nij)*(b(j)-nij)/(nij+1)/(N-a(i)-b(j)+nij+1);  
        end
        CC=N*(a(i)-1)*(b(j)-1)/a(i)/b(j)/(N-1)+N/a(i)/b(j);
        bound(i,j)=a(i)*b(j)/N^2*log2(CC);         
    end
end

EMI_bound=sum(sum(bound));
EMI_bound_2=log2(r*c/N+(N-r)*(N-c)/(N*(N-1)));
EMI=sum(sum(EPLNP-E3));

AMI_=(MI-EMI)/(sqrt(Ha*Hb)-EMI);
NMI=MI/sqrt(Ha*Hb);


%If expected mutual information negligible, use NMI.
if abs(EMI)>EMI_bound
    fprintf('The EMI is small: EMI < %f, setting AMI=NMI',EMI_bound);
    AMI_=NMI;
end;

%---------------------auxiliary functions---------------------
function Cont=Contingency(Mem1,Mem2)

if nargin < 2 || min(size(Mem1)) > 1 || min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end