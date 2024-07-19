function [pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,del,forFeedback,coeff,beProps,alphaIn)%%**
%% Behaviour: default "be" created but not used (eval calculation has expandB=0)
%% Kroenecker products of all necesary parameters and initial D done at the end - still to do

%Fix vax rollout
data.atimes=[343,367,412];
data.arates(3)=3*data.arates(3);
%data.arates=[];

%plotSingleRun(dataOcc,ones(1,19),dataUK1,[poptimGVAstart,.2],Xin20ld,coeff(:,1:3)');
%% MODEL TYPE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
feedbackIn=1;
numPCA=size(coeff,1);%Number of principle components - still manually specify parameters under "feedback"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%data.NNs - column vector of population
%Possible generalsiation to within-sector heterogeneity - one column per subsector

%% POPULATION PARAMETERS:

%Population Density
[n,na]=size(data.NNs);
ntot=n*na;
NN=sum(data.NNs,2);
NNbar=reshape(data.NNs,ntot,1);
NNrep=repmat(NN,na,1);

%Age-Sector Breakdown
lc=4;
adInd=3;
lx=length(data.NNs)-lc;
xin=data.Xdata(1:lx).^(1./data.alpha);%*Phew
xin=xin';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Default (no behavioural change) settings:
be=struct;
be.alphaB=1;
be.propsB=beProps;%0*ones(5,1);%Column vector of proportions in each age group
be.expandB=0;%For eigenvalue calculation
be.numGroups=2;%Hard-coded in "makeDs" for now
be.nonWorkInds=(lx+2)*be.numGroups+(1:be.numGroups);%lx overwritten later to include be groups
be.NNfull=kron(NNbar,ones(2,1));%Age/sector population, irrespective of behavioural group

%D=beMakeDs(NN,xin,data,zeros(1,lx),be,zeros(5,1));%propsB=zeros(5,1)
D=beMakeDs(NN,1,data,zeros(1,lx),be,zeros(5,1));%propsB=zeros(5,1)

Dout=D;%Now redundant
be.alphaB=alphaIn(1);%1;
be.alphaB2=alphaIn(2);

be.alphaB3=alphaIn(3);

be.expandB=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PCA feedback:

pr=struct;
%Last of each distinct x:
%pr.x20=[3     4    23    24    26    27    30    35    36    40    43    44    49    53    54    55    57    59    62   63]';
pr.x20=1:size(coeff,2);%Factor analysis/22 sectors in
%
%Number of each x %Factor analysis:
pr.xrep=[3     1    19     1     2     1     3     5     1     4     3     1     5     4     1     1     2     2     3     1]';
pr.xrep=ones(1,20);%Factor analysis/22 sectors in

%% Coefficients for dimensional reduction
pr.coeff=coeff;
%pr.coeff=zeros(4,22); pr.coeff(1,11)=1; pr.coeff(2,18)=1; pr.coeff(3,20)=1; pr.coeff(4,10)=1;%numPCA sectors isolated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pr.a=data.alpha;%0.411475;
%Track and trace:
pr.tntOn=426;
pr.redSus=0;
pr.redAsy=0;
pr.redSym=0;
pr.logSectors=62;
pr.optimFrom=1;

%Feedback parameters:
if feedbackIn==1
    %pr.coeff=coeff(:,1:numPCA)';%PCA
    %pr.coeff=coeff;%Factor analysis
    %
    pr.L1=forFeedback(1);%0.4;
    %
    pr.k1=forFeedback(2:2+numPCA-1)';%(2);%.0005;
    pr.H01=forFeedback(numPCA+2);%(3);%1e4;
    %pr.m1=forFeedback(numPCA+3);
    %}
    %pr.betamatrix=reshape(forFeedback,5,3);
    %}
else
    pr.L1=0;%0.4;
    pr.k1=1;%.0005;
    pr.H01=10e5;%1e4;
    pr.m1=1;
end

%% DISEASE PARAMETERS:

%from Knock et al. (2020) unless indicated otherwise

%Transmission
pr.R0=R0;%2.6882;%from fitting
pr.red=0.58;%from Byambasuren et al. (2020)

%Latency and Onset
Text=4.6;
pr.sigma=1/Text;
%Tonset=1;
%pr.omega=1/Tonset;

%Case Pathways
pr.p1=0.6;
[ph,pd,~,~]=heParamsAge(data,pr.p1);
ph=[repmat(ph(adInd),lx,1);ph];
pd=[repmat(pd(adInd),lx,1);pd];

%Hospitalisation and Death Rates (proportion and rate combined)
Tsh=4.5;%from Global
Thd=data.Thd;%14; %DH - as Sri Lanka

%Recovery Rates (proportion and rate combined)
Ta=2.1;%asymptomatic
pr.g1=1/Ta;
Ts=4;%symptomatic
Threc=data.Threc;%hospitalised%14;

%DH:
time_symp=ph.*Tsh+(1-ph).*Ts;
pr.g2=(1-ph)./time_symp;
pr.h=ph./time_symp;
time_hosp=pd.*Thd+(1-pd).*Threc;
pr.g3=(1-pd)./time_hosp;
pr.mu=pd./time_hosp;
%Delta variant:
pr.hDelta=1.85*pr.h;
pr.g2Delta=(1-1.85*ph)./time_symp;
pr.g3Delta=(1-pd)/7.6;
pr.muDelta=pd/7.6;

%Immunity Loss
Ti=365;%from Global
pr.nu=1/Ti;

%% BETA CALIBRATION:

Deff=Dout.*repmat(NNbar,1,n*na)./repmat(NNbar',n*na,1);
onesn=ones(ntot,1);

F=zeros(3*ntot,3*ntot);
F(1:ntot,ntot+1:end)=[pr.red*Deff,Deff];

vvec=[pr.sigma.*onesn;      pr.g1.*onesn;       (pr.g2+pr.h).*onesn];%g2 and h are vectors
V=diag(vvec);
V(ntot+1:2*ntot,1:ntot)=    diag(-(1-pr.p1) .*pr.sigma  .*onesn);
V(2*ntot+1:3*ntot,1:ntot)=  diag(-pr.p1     .*pr.sigma  .*onesn);

GD=F/V;

d=eigs(GD,1);%largest in magnitude (+/-) 
R0a=max(d); 
beta=pr.R0/R0a;%beta scales actual R0 to fitted R0
pr.betaAddOn=beta*(5.08/2.8025-1);

%% PREPAREDNESS PARAMETERS:

%Adherence to NPIs
pr.betamod=[1,1,del];

%Self-Isolation and Quarantine
pr.p3=0.00;%proportion of asymptomatic self-isolating
pr.p4=0.00;%proportion of symptomatic self-isolating
%pr.p5=0.00;%proportion of severe self-isolating
%
pr.odds=0;%asymptomatic quarantining rate
pr.q1=0;%mild quarantining rate
pr.q2=0;%severe quarantining rate


%Hospital Capacity
pr.Hmax=data.Hmax;%2017;

pdoc=heParamsAge_oc(data,pr.p1);
pdoc=[repmat(pdoc(adInd),lx,1);pdoc];
pr.mu_oc=pdoc/Thd;
pr.g3_oc=(1-pdoc)/Threc;

%% VACCINATION PARAMETERS:

vx=struct;%from RTM unless indicated otherwise

%Vaccine 1 (two doses)
vx.hrv1=    1/28;                       %time to develop v-acquired immunity (AstraZeneca)
vx.scv1=    0.58;                       %infection-blocking efficacy
vx.scv1b=   0.4;                        %delta variant infection-blocking efficacy
vx.p1v1=    0;                          %disease-blocking efficacy          
vx.hv1=     1-((1-0.90)/(1-vx.scv1));   %severe-disease-blocking efficacy
vx.dv1=     0;                          %death-blocking efficacy
vx.trv1=    0.40;                       %transmission-blocking efficacy
vx.nuv1=    1/Inf;                      %duration of v-acquired immunity

vx.h_v1=    (1-vx.hv1)*ph/Tsh;
vx.g2_v1=   (1-(1-vx.hv1)*ph)/Ts;
vx.mu_v1=   (1-vx.dv1)*pd/Thd;
vx.g3_v1=   (1-(1-vx.dv1)*pd)/Threc;

% vx.mu_ocv1= 1.19*vx.mu_v1;%Wilde et al. (2021)
% vx.g3_ocv1= (1-1.19*(1-vx.dv1)*pd)/Threc;
vx.mu_ocv1=   (1-vx.dv1)*pdoc/Thd;
vx.g3_ocv1=   (1-(1-vx.dv1)*pdoc)/Threc;

%Administration Rates
ta=data.atimes;%[391,523,548];
arate=1*data.arates;%[9200,16666,112254];

%Population Uptake
Npop=   data.Npop;
NNage=  Npop';%DH %[Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:end))];
pr.up=  data.uptake;%[0.00,0.52*0.90,0.90,0.90];%12:0.52, 14:0.39, 16:0.26, 18:0.13
NNwv=   NNage.*pr.up;

%% Correction to vax rates:
%
vx.startp1=345;
vx.startp2=411;
vx.startp3=579;
vx.startp4=580;
vx.startp5=581;
vx.end=582;

vx.aratep1=[0;0;0;0;170050];
%vx.aratep2=[0;0;2.1296e+05;0];
vx.aratep2=[1.6952e+05;0;0;4.3438e+04;0];
vx.aratep3=zeros(5,1);
vx.aratep4=zeros(5,1);
vx.aratep5=zeros(5,1);
%}
%% Kron for behaviour:
kn=ones(2,1);
%k4=ones(2,1);
pr.g2=kron(pr.g2,kn);
pr.h=kron(pr.h,kn);
pr.g3=kron(pr.g3,kn);
pr.mu=kron(pr.mu,kn);
pr.mu_oc=kron(pr.mu_oc,kn);
pr.g3_oc=kron(pr.g3_oc,kn);

pr.g2Delta=kron(pr.g2Delta,kn);
pr.hDelta=kron(pr.hDelta,kn);
pr.g3Delta=kron(pr.g3Delta,kn);
pr.muDelta=kron(pr.muDelta,kn);

vx.h_v1=kron(vx.h_v1,kn);
vx.g2_v1=kron(vx.g2_v1,kn);
vx.mu_v1=kron(vx.mu_v1,kn);
vx.g3_v1=kron(vx.g3_v1,kn);
vx.mu_ocv1=kron(vx.mu_ocv1,kn);
vx.g3_ocv1=kron(vx.g3_ocv1,kn);
%4 long:
vx.aratep1=kron(vx.aratep1,kn);
vx.aratep2=kron(vx.aratep2,kn);
vx.aratep3=kron(vx.aratep3,kn);
vx.aratep4=kron(vx.aratep4,kn);
vx.aratep5=kron(vx.aratep5,kn);
end

function [phgs,pdgh,Threc,Thd]=heParamsAge(datax,ps)
%%
nn=datax.Npop';

ranges=[1,3,9,4];
%nntot=repelem(nn,ranges);%DH
nnprop=repelem(1./ranges,ranges);%DH

subs=1:4;
subs=repelem(subs,ranges);

%%

ihr=    [0.030000 0.002600 0.000840	0.000420 0.000800 0.002600 0.004000	0.006300 0.012000 0.019000 0.023000	0.040000 0.096000 0.100000 0.240000	0.500000 0.500000];
phgs=   ihr./ps;

ifr=    [0.000310 0.000030 0.000010 0.000000 0.000010 0.000040 0.000060 0.000130 0.000310 0.000700 0.001160 0.002760 0.008670 0.012150 0.035120 0.084300 0.096960];
pdgh=   ifr./ihr;

phgs=   accumarray(subs',phgs.*nnprop);
pdgh=   accumarray(subs',pdgh.*nnprop);

Threc=  NaN;
Thd=    NaN;

end

function pdgh=heParamsAge_oc(datax,ps)

nn=datax.Npop';
nn=datax.Npop';

ranges=[1,3,9,4];
%nntot=repelem(nn,ranges);%DH
nnprop=repelem(1./ranges,ranges);%DH

subs=1:4;
subs=repelem(subs,ranges);

ihr=    [0.030000 0.002600 0.000840	0.000420 0.000800 0.002600 0.004000	0.006300 0.012000 0.019000 0.023000	0.040000 0.096000 0.100000 0.240000	0.500000 0.500000];

ifr=    [0.0003 0.0003 0.0002 0.0002 0.0003 0.0003 0.0009 0.0009 0.0034 0.0034 0.0073 0.0073 0.0253 0.0253 0.064 0.064 0.1325];
pdgh=   ifr./ihr;

pdgh=   accumarray(subs',pdgh.*nnprop);

end