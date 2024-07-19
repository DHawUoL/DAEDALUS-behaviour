function f=heMakeDs(NN,x,data,wfh,be,propsB)
%% Behaviour: end section only. 2 classes hard coded in this script

%NN - vector of populations including non-working.
%x - proportion of each sector open - not including non-working.
%The contact matrices are entirely determined by 'economics', both sector closures and working-from-home
%However, the succeptability of children may be included

%% COMMUNITY-COMMUNITY MATRIX:

childrenLessSus=0;

C=data.CM;%source: polymod-uk-all

%%

na=length(C);
lx=length(x);%Number of sectors%make this standard
ln=length(NN);

adInd=3;%Adult index
CworkRow=C(adInd,:);

NNrep=repmat(NN'/sum(NN),ln,1);%population proportion matrix
NNrel=NN([1:lx,lx+adInd])/sum(NN([1:lx,lx+adInd]));%proportion of adult population in each sector as vector

%Make A:
matA=zeros(ln,ln);
matA(lx+1:end,lx+1:end)=C;
matA(1:lx,lx+1:end)=repmat(CworkRow,lx,1);
matA(:,[1:lx,lx+adInd])=repmat(matA(:,lx+adInd),1,lx+1).*repmat(NNrel',ln,1);

%% Modify depending on x:

if lx==1
    %
    %Education:
    matA(lx+1,lx+1)=    matA(lx+1,lx+1)+    x^2*  data.schoolA1;%mixing within age groups only
    matA(lx+2,lx+2)=    matA(lx+2,lx+2)+    x^2*  data.schoolA2;
    
    %Hospitality:
    %sects=[36,58,59,60,62];
    %psub=data.NNs(sects);
    psub=data.fracHosp*x^2;%sum(psub.*x(sects).^2)/sum(psub);%constant from 0-1, weighted measure of how much sectors are open
    matA([1:lx,lx+adInd],:)=    matA([1:lx,lx+adInd],:)+    psub*   data.hospA3*    NNrep([1:lx,lx+adInd],:);%mixing between all age groups, including pre-school
    matA(lx+2,:)=               matA(lx+2,:)+               psub*   data.hospA2*    NNrep(lx+2,:);
    matA(ln,:)=                 matA(ln,:)+                 psub*   data.hospA4*    NNrep(ln,:);
    %}
elseif lx==10
    %Education:
    matA(lx+1,lx+1)=matA(lx+1,lx+1)+data.propschools*data.schoolA1*x(9);
    matA(lx+2,lx+2)=matA(lx+2,lx+2)+data.propschools*data.schoolA2*x(9);
    
    %Hospitality:
    %matA([1:lx,lx+3:ln],:)=matA([1:lx,lx+3:ln],:)+NNrep([1:lx,lx+3:ln],:)*datax.prophosp*x(10)*datax.hospA34(1);
    matA([1:lx,lx+3],:)=matA([1:lx,lx+3],:)+NNrep([1:lx,lx+3],:)*data.prophosp*x(10)*data.hospA3(1);
    matA(ln,:)=matA(ln,:)+NNrep(ln,:)*data.prophosp*x(10)*data.hospA4(1);
    matA(lx+2,:)=matA(lx+2,:)+NNrep(lx+2,:)*data.prophosp*x(10)*data.hospA2;%(1)
      
elseif lx==35
    %Education:
    matA(lx+1,lx+1)=    matA(lx+1,lx+1)+    x(32)*  data.schoolA1;%mixing within age groups only
    matA(lx+2,lx+2)=    matA(lx+2,lx+2)+    x(32)*  data.schoolA2;
    
    %Hospitality:
    sects=[22,34];
    psub=data.NNs(sects);
    psub=sum(psub.*x(sects))/sum(psub);%constant from 0-1, weighted measure of how much sectors are open
    matA([1:lx,lx+adInd],:)=    matA([1:lx,lx+adInd],:)+    psub*   data.hospA3*    NNrep([1:lx,lx+adInd],:);%mixing between all age groups, including pre-school
    matA(lx+2,:)=               matA(lx+2,:)+               psub*   data.hospA2*    NNrep(lx+2,:);
    matA(ln,:)=                 matA(ln,:)+                 psub*   data.hospA4*    NNrep(ln,:);

elseif lx==36
    %Education:
    matA(lx+1,lx+1)=    matA(lx+1,lx+1)+    x(33)*  data.schoolA1;%mixing within age groups only
    matA(lx+2,lx+2)=    matA(lx+2,lx+2)+    x(33)*  data.schoolA2;
    
    %Hospitality:
    sects=[25,35];
    psub=data.NNs(sects);
    psub=sum(psub.*x(sects))/sum(psub);%constant from 0-1, weighted measure of how much sectors are open
    matA([1:lx,lx+adInd],:)=    matA([1:lx,lx+adInd],:)+    psub*   data.hospA3*    NNrep([1:lx,lx+adInd],:);%mixing between all age groups, including pre-school
    matA(lx+2,:)=               matA(lx+2,:)+               psub*   data.hospA2*    NNrep(lx+2,:);
    matA(ln,:)=                 matA(ln,:)+                 psub*   data.hospA4*    NNrep(ln,:);

elseif lx==63
    %Education:
    matA(lx+1,lx+1)=    matA(lx+1,lx+1)+    x(55)^2*  data.schoolA1;%mixing within age groups only
    matA(lx+2,lx+2)=    matA(lx+2,lx+2)+    x(55)^2*  data.schoolA2;
    
    %Hospitality:
    sects=[36,58,59,60,62];
    psub=data.NNs(sects);
    psub=sum(psub.*x(sects).^2)/sum(psub);%constant from 0-1, weighted measure of how much sectors are open
    matA([1:lx,lx+adInd],:)=    matA([1:lx,lx+adInd],:)+    psub*   data.hospA3*    NNrep([1:lx,lx+adInd],:);%mixing between all age groups, including pre-school
    matA(lx+2,:)=               matA(lx+2,:)+               psub*   data.hospA2*    NNrep(lx+2,:);
    matA(ln,:)=                 matA(ln,:)+                 psub*   data.hospA4*    NNrep(ln,:);

else
    error('Unknown economic configuration!');
    
end

%Transport:
matA(1:lx,1:lx)=    matA(1:lx,1:lx)+    repmat(x',lx,1).*   data.travelA3.*  NNrep(1:lx,1:lx).*  repmat(1-wfh,lx,1).*repmat(1-wfh',1,lx);%home-working has a compound effect

%% WORKER-WORKER AND COMMUNITY-WORKER MATRICES:

%Make B and C:
valB=data.B;
valB=valB.*(1-wfh).*(1-wfh);%home-working has a compound effect
valC=data.C;
valC=valC.*(1-wfh);
valB(lx+1:ln)=0;
valC(lx:1:ln)=0;

%Note: bug here in single sector case (at least) - fixed for behaviour
xx=[x;zeros(ln-lx,1)];%x(lx+1:ln,1)=0; %*1sector

matB=diag(xx.*valB');
matC=repmat(xx.*valC',1,ln).*NNrep;

%%p

D=matA+matB+matC;%dot(sum(matB+matC,2),data.NNs)/sum(data.NNs)

%%

if childrenLessSus==1
    D(end-3,:)=(1-0.5)*D(end-3,:);
    D(end-2,:)=(1-(0.5*11/15))*D(end-2,:);
end

%f=D;

%% Expand D out by 2 behaviour class - BE
% - assumes 2 behaviour groups and can further split the working age class
% by sector
if be.expandB==1
    numD=length(D);
    %behMat=repmat([be.alphaB^2,be.alphaB;be.alphaB,1],numD,numD);
    behMat=[be.alphaB^2,be.alphaB;be.alphaB,1];
    pOrder=reshape([propsB';1-propsB'],2*length(propsB),1);%be.probsB has length number of age/sector compartments
    pOrder=repmat(pOrder',numD*2,1);
    
    %behmat=repmat([propsB',1-propsB'],numD,1);%pr.propsB - column vector of proportions in each age group
    %behmat=repmat(behmat,length(D),1);
    %D=kron(D,ones(2)).*behMat.*pOrder;
    D=kron(D,behMat).*pOrder;
end
f=D;
%{
Parameters to add:
pr.alphaB
pr.propsB - column vector of proportions in each age group
%}
end

function out = inverse_kron(K,array,input_id)

switch input_id
    case 1
        out = K(1:size(K,1)/size(array,1),1:size(K,2)/size(array,2))./array(1);
    case 2
        out = K(1:size(array,1):end,1:size(array,2):end)./array(1);
    otherwise
        error('The Input ID must be either 1 or 2')
end
end