function [f,g,tvec,rhohat]=beSimCovid19vax(pr,be,vx,beta,tvec,Xit,ntot,NNvec,phi1,phi2,seedvec,S0,tau,plotTau,pOrder,data)

%% PARAMETERS/MODEL TYPE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
feedbackIn=0;%=0 for variant with behavioural groups
lx=be.numGroups;%%BH
NNbar=kron(NNvec(:,1),[0;1]);%%BH - NNvec stays aggregated by behavioural group
S0=NNbar;
%{
numSect=lx/be.numGroups;
NNbar=NNvec(:,1);
sumWorkingAge=sum(NNbar([1:lx,lx+(2*be.numGroups+1:3*be.numGroups)]));%sum(NNbar([1:lx,lx+3])); %PD check
%}
%%BH:
NNvecHold=repmat(zeros(size(NNvec)),2,1);
%NNvecHold=repmat(NNvecHold,2,1);

nc=20;%Number of model compartments
solvetype=2;
hospInc=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For behaviour - assumes 2 behaviour groups
indsForTotals=kron((2:2:ntot)',[1;0])+kron((1:2:ntot-1)',[0;1]);

%% DETERMINISTIC OR STOCHASTIC SOLVER:

if solvetype==2
    
    lt=length(tvec);
    zn=zeros(ntot,1);
    
    t0=tvec(1);
    y0=[S0;repmat(zn,6,1);NNbar-S0;repmat(zn,nc-8,1)];
    
    toutAll=t0;
    Sout=S0';
    Sv1out=repmat(zn',1,2);
    %Iout=sum(seedvec');
    Iasy=zn';
    Isym=zn';%Should be seed - but will only affect first time point
    Iasyvax=zn';
    Isymvax=zn';
    Hout=zn';
    HnewAll=[];
    Dout=zn';
    DEout=repmat(zn,1,lt);
    %Rout=sum((NNbar-S0)');
    Vout=zn';
    incOut=zn';
    hospIncOut=0;
    Rt=zeros(lt-1,1);
    
    %% FIXED INPUT:
    
    Houtlast=0;
    pr.rhohat=zeros(1,lt);
    pr.Hout=zeros(1,lt-1);
    for i=1:lt-1
        
        %% MODIFIERS A:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%**
        %betaIn=beta;%Single fitted modifier
        %pr.rhohat(i)=1;
        %
        if feedbackIn==1
            %
            if i<3
                betaIn=beta*ones(5,1);%Single fitted modifier %Beta vector
                pr.rhohat(i)=1;
            else
                pr.rhohat(i)=rho2ofh(pr.coeff*(pr.xfull(pr.x20,i)-pr.Xmu),pr);
                betaIn=beta*pr.rhohat(i);
            end
        else
            %*pr.betamod(i) is included in making of D's
            betaIn=beta;%Single fitted modifier
            pr.rhohat(i)=1;
            %}
            %%BH:
            if i<4%******** 2 for initial fit, <4 for full fit
                propsBi=zeros(5,1);
            else
                propsBi=be.propsB(:,i);%BH be.propsB(:,i)=propsBi;
            end
                %}

            if tvec(i)>=183%183%i>32%17%******** 1 for initial fit, 4 for full fit 14
                be.alphaB=be.alphaB2;%Overwrite - used in Dvec=
            end
        
            if i>100%12
                be.alphaB=be.alphaB3;%Overwrite - used in **
            end

            %%BH Dvec(:,:,i)=pr.betamod(i)*beMakeDs(NNvec(:,i),Xit(:,i),data,data.wfhAv(i,:),be,propsBi);%be.propsB(:,i));
            pOrderi=reshape([propsBi';1-propsBi'],ntot,1);
            NNvecHold(:,i)=kron(NNvec(:,i),ones(2,1)).*pOrderi;
            pOrder(:,i)=pOrderi;
            D=pr.betamod(i)*beMakeDs(NNvec(:,i),Xit(:,i),data,data.wfhAv(i,:),be,propsBi);
            
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %
        %Delta variant
        %
        if tvec(i)>=487 && tvec(i-1)<487 %i==17
            %Overwrite once:
            beta=beta+.5*pr.betaAddOn;%del=1 always
            pr.h=pr.hDelta;
            pr.g2=pr.g2Delta;
            pr.g3=pr.g3Delta;
            pr.mu=pr.muDelta;
        elseif tvec(i)>487 %i==18
            beta=beta+.5*pr.betaAddOn;
        end
        %}
        tend=tvec(i+1);
        NNfeed=NNvecHold(:,i);%%BH
        NNfeed(NNfeed==0)=1;%NNfeed functions as NN0, i.e. zeros changed to ones
        %%BH D=Dvec(:,:,i);
        
        %Vaccination Rollout by Sector
        %
        vx.ratep1period=vx.aratep1.*pOrder(:,i);%Does all vax periods - change to loop;
        vx.ratep2period=vx.aratep2.*pOrder(:,i);
        vx.ratep3period=vx.aratep3.*pOrder(:,i);
        vx.ratep4period=vx.aratep4.*pOrder(:,i);
        vx.ratep5period=vx.aratep5.*pOrder(:,i);
        %}
        %% DH commented out for hospital admissions fit
        %{
        if hospInc>=1
            Rt(i)=heComputeEigs(pr,beta,D,NNfeed,ntot,y0(1:lx+4));
        end
        %}
        
        seedvecIn=seedvec.*pOrder(:,i);
        [tout,Sclass,Hclass,Dclass,DEcum,Rcum,y0,Hnew,Sv1class,Vclass,Iasyclass,Isymclass,Iasyvaxclass,Isymvaxclass,incOutclass,hospIncOutclass]=integr8(pr,vx,betaIn,ntot,NNfeed,NNvecHold(:,i),D,phi1,phi2,seedvecIn,t0,tend,y0,i,hospInc);
        
        toutAll=[toutAll;tout(2:end)];
        Sout=[Sout;Sclass(2:end,:)];
        Sv1out=[Sv1out;Sv1class(2:end,:)];
        %Iout=[Iout;Itot(2:end)];
        Hout=[Hout;Hclass(2:end,:)];
        HnewAll=[HnewAll;Hnew(2:end)];
        Dout=[Dout;Dclass(2:end,:)];
        DEout(:,i)=DEcum;
        %Rout(:,i)=Rcum;
        Vout=[Vout;Vclass(2:end,:)];
        Iasy=[Iasy;Iasyclass(2:end,:)];
        Isym=[Isym;Isymclass(2:end,:)];
        Iasyvax=[Iasyvax;Iasyvaxclass(2:end,:)];
        Isymvax=[Isymvax;Isymvaxclass(2:end,:)];
        incOut=[incOut;incOutclass(2:end,:)];
        %Not part of y:
        hospIncOut=[hospIncOut;sum(hospIncOutclass(2:end,:),2)];

        Houtlast=mean(sum(Hclass,2));
        pr.Hout(i)=Houtlast;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Change behavioural parameters based on occupancy:
        %{
        OccNow=sum(Hout(end,:));

        if i<3
            propsBi=zeros(5,1);
        else
            propsBi=beFeedback(pr.coeff*(pr.xfull(pr.x20,i)-pr.Xmu),pr)*ones(5,1);
            %propsBi=be.propsB(:,i);
        end
        be.propsB(:,i)=propsBi;
        
        Dvec(:,:,i)=pr.betamod(i)*beMakeDs(NNvec(:,i),Xit(:,i),data,data.wfhAv(i,:),be,propsBi);%be.propsB(:,i));
        pOrderi=reshape([be.propsB(:,i)';1-be.propsB(:,i)'],ntot,1);
        NNvecHold(:,i)=kron(NNvec(:,i),ones(2,1)).*pOrderi;
        pOrder(:,i)=pOrderi;
        NNvec=NNvecHold;%%xx Not ncessarily integer populations
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        t0=tend;
        if i<lt-1
            %% Shift behavioural classes first:
            % - assumes 2 behaviour groups
            
            %Fitting link function:
            %


            if i<3%******** 3 for main model 2 for initial fit
                propsBi=zeros(5,1);
            else %Already have i<lt-1
                %
                Hx1=sum(Hout(end,:));
                Hx2=sum(HnewAll(end,:));
                Hx3=Hx1-sum(Hout(end-13,:));
                %propsBi=beFeedback2([Hx1,Hx2,Hx3],pr)*ones(5,1);%ones(5,1);%Feedback function
                propsBi=beFeedback2([pr.xfull(:,i+1)',Hx1/1e5],pr)*ones(5,1);%ones(5,1);
                %propsBi=beFeedback2([pr.xfull(:,i+1)'],pr)*ones(5,1);%ones(5,1);
                %}
                %propsBi=be.BiFirstFit*ones(5,1);%******** Line for initial fit
            end

            
            %}
            %Fitting individual p's:
            %{
            if i<3%********
                propsBi=zeros(5,1);
            else
                propsBi=pr.xfull(i+1)*ones(5,1);
            end
            %}
            be.propsB(:,i+1)=propsBi;

            beDiff=be.propsB(:,i+1)-be.propsB(:,i);%From risk to risk averse
            beDiff=[beDiff';-beDiff'];
            %beDiffOut=-beDiff;
            %beDiff(beDiff<0)=0;
            beDiff=reshape(beDiff,ntot,1);
            beDiff=repmat(beDiff,1,nc);
            %From these behavioural classes:
            %beDiffOut(beDiffOut<0)=0;
            %beDiffOut=reshape(beDiff,ntot,1);
            %beDiffOut=repmat(bdiffOut,1,nc);
            %
            y0=reshape(y0,[ntot,nc]);%IC
            y0totals=y0+y0(indsForTotals,:);%Totals in each behavioural group
            y0=y0+beDiff.*y0totals;
            y0=reshape(y0,ntot*nc,1);
        end
    end
    pr.rhohat(end)=yofh(Houtlast);
    Iout=sum(Iasy+Isym,2);
    incOut=diff(diff(incOut(end-2:end)));
    HdiffOut=diff(Hout(end-1:end));
    
    %% SWITCHING in here
    
elseif solvetype==1
    
    error('Final size calculations not possible')
    
elseif solvetype==3
    
    error('Code not written yet')
    %f=stochSim(y0,beta,gamma,n,ntot,NN,NN0,D,seed,phi1,phi2,tau,alpha);
    
end

%% OUTPUTS:
 
if tau==plotTau

    dodiff=0;%Can't set to 1 if have waning immunity
    plotEpi(toutAll,Iout,Hout,ntot,dodiff,tvec);%-tvec(2)+38

end    

if hospInc==0
    
    f=[toutAll,...
       sum(Sout,2),...
       Iout,...
       sum(Hout,2),...
       sum(Dout,2),...
       sum(Vout(:,lx+1:lx+be.numGroups),2),...
       sum(Vout(:,lx+be.numGroups+1:lx+2*be.numGroups),2),...
       sum(Vout(:,[1:lx,lx+2*be.numGroups+1:lx+3*be.numGroups]),2),...
       sum(Vout(:,lx+3*be.numGroups+1:lx+4*be.numGroups),2)];%sum(DEout(end,:));
    %f=[toutAll(toutAll>0),sum(Hout(toutAll>0,:),2)];%For epi fit
    %Main constraints
    f(toutAll<1,:)=[];
    %For track and trace:
    %g=[Iasy,Isym,Iasyvax,Isymvax];
    %g(toutAll<1,:)=[];
    %For optimisation:
    g=[max(f(pr.optimFrom:end,4)),...
       HdiffOut];%Rt(end)];
   %Rt calculated at end of each period
   rhohat=pr.rhohat;%pr.Hout; pr.rhohat;
   
elseif hospInc==1%****
    
    f=[toutAll,...
       sum(Sout,2),...
       Iout,...
       sum(Hout,2),...
       sum(Dout,2),...
       sum(Vout(:,lx+1),2),...
       sum(Vout(:,lx+2),2),...
       sum(Vout(:,[1:lx,lx+3]),2),...
       sum(Vout(:,lx+4),2)];
    f(toutAll<1,:)=[];

    %g=[max(f(pr.optimFrom:end,4)),...
    %   HdiffOut];%Rt(end)];
    g=hospIncOut;
    g(toutAll<1,:)=[];

   rhohat=be.propsB(1,:);%pr.rhohat;%pr.Hout; pr.rhohat;
    
elseif hospInc==2
    
    f=[toutAll(toutAll>0),sum(Sout(toutAll>0,:),2),sum(Hout(toutAll>0,:),2)];%Occupancy
    g=sum(sum(HnewAll(toutAll>tvec(3),:)));%Incidence
    rhohat=pr.rhohat;
    
end

end

function f=beFeedback2(H,pr)
f=pr.L1/(1+exp(dot(pr.k1,H)-pr.H01'));
%f=pr.L1*normcdf(dot(pr.k1,H)-pr.H01');
end

%%

function [tout,Sclass,Hclass,Dclass,DEcum,Rcum,y0new,Hnew,Sv1class,Vclass,Iasyclass,Isymclass,Iasyvaxclass,Isymvaxclass,incidence,hospIncOut]=integr8(pr,vx,beta,ntot,NN0,NNbar,D,phi1,phi2,seedvec,t0,tend,y0,i,hospInc)
%% CALL:

fun=@(t,y)integr8covid(t,y,pr,vx,beta,ntot,NN0,NNbar,D,phi1,phi2,seedvec,i);
%{
if i>1 && t0~=ceil(t0)
    [tout,yout]=ode45(fun,[t0,ceil(t0),ceil(t0)+1:floor(tend)],y0);
    tout(1)=[];
    yout(1,:)=[];
else
    [tout,yout]=ode45(fun,[t0,ceil(t0)+1:floor(tend)],y0);
end
%}
[tout,yout]=ode45(fun,[t0,tend],y0);
    
     Hnew=0;

%% EC:

Sclass=     yout(:,  0*ntot+1:1*ntot);
Sv1class=   yout(:,  8*ntot+1:10*ntot);

Iasyclass=  yout(:, 2*ntot+1:3*ntot)    +yout(:, 3*ntot+1:4*ntot);
Iasyvaxclass=  yout(:, 11*ntot+1:12*ntot)    +yout(:, 12*ntot+1:13*ntot);
Isymclass=  yout(:, 4*ntot+1:5*ntot)    +yout(:, 5*ntot+1:6*ntot);
Isymvaxclass=  yout(:, 13*ntot+1:14*ntot)    +yout(:, 14*ntot+1:15*ntot);
%Itot=   sum(yout(:,  2*ntot+1:6*ntot),2)+sum(yout(:,  11*ntot+1:15*ntot),2);

Hclass=     yout(:,  6*ntot+1:7*ntot)       +yout(:,  15*ntot+1:16*ntot);
Rcum=       yout(end,7*ntot+1:8*ntot)       +yout(end,16*ntot+1:17*ntot);
Dclass=     yout(:,  17*ntot+1:18*ntot);
DEcum=      yout(end,17*ntot+1:18*ntot);
Vclass=     yout(:,18*ntot+1:19*ntot);
incidence=  yout(:,19*ntot+1:20*ntot);
y0new=      yout(end,:)';

if hospInc>=1
    %[~,hospIncOut]=fun(tout,yout);
    yy=yout;%(tout>=0,:);
    tt=tout;%(tout>=0);
    [~,hospIncOut] = cellfun(fun, num2cell(tt), num2cell(yy',1)','uni',0);
    hospIncOut=cell2mat(hospIncOut);
else
    hospIncOut=0;
end

end

%%

function [f,g]=integr8covid(t,y,pr,vx,betaIn,ntot,NN0,NNbar,D,phi1,phi2,seedvec,period)
%% IC:

S=      y(0*ntot+1:1*ntot);
E=      y(1*ntot+1:2*ntot);
Ina=    y(2*ntot+1:3*ntot);
Isa=    y(3*ntot+1:4*ntot);
Ins=    y(4*ntot+1:5*ntot);
Iss=    y(5*ntot+1:6*ntot);
H=      y(6*ntot+1:7*ntot);
R=      y(7*ntot+1:8*ntot);

Shv1=   y(8*ntot+1:9*ntot);
Sv1=    y(9*ntot+1:10*ntot);
Ev1=    y(10*ntot+1:11*ntot);
Inav1=  y(11*ntot+1:12*ntot);
Isav1=  y(12*ntot+1:13*ntot);
Insv1=  y(13*ntot+1:14*ntot);
Issv1=  y(14*ntot+1:15*ntot);
Hv1=    y(15*ntot+1:16*ntot);
Rv1=    y(16*ntot+1:17*ntot);

DE=     y(17*ntot+1:18*ntot);
V=      y(18*ntot+1:19*ntot);

%% FOI:
%% MODIFIERS B:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%phi=phi1;
phi=phi1+phi2*cos(pi*(t-pr.tlag)/180);%Seasonality****

%%**
beta=betaIn;
%{
if period<3
    beta=betaIn;
else
    %beta=betaIn*rho2ofh(sum(H),pr);
    %beta=betaIn*rho2ofh([pr.xfull([53,43,49],period);sum(H)],pr);%[36,58,60]
    %beta=betaIn*rho2ofh([pr.coeff*pr.xfull(pr.x20,period);sum(H)],pr);

    Hin=pr.xfull*pr.coeff; Hin=Hin(period,:)';%Factor analysis
    beta=betaIn*rho2ofh([Hin;sum(H)],pr);%Factor analysis
    %beta=betaIn*rho2ofh(pr.coeff*pr.xfull(pr.x20,period),pr);
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Track and trace:
I=(pr.red*Ina+Ins)  +(1-vx.trv1)*(pr.red*Inav1+Insv1);%Only non-self-isolating compartments
redSus=1;

foi=phi*beta.*(D*(I./NN0));

if period<3
    seed=phi*(seedvec./NN0);
else
    seed=zeros(ntot,1);
end

%% HOSPITAL OCCUPANCY:

occ=max(1,sum(H+Hv1));

%% SELF-ISOLATION:
    
p3=0;
p4=0;

%% VACCINATION:

nonVax=NNbar-V;
nonVax(nonVax==0)=1;
%S (and R) ./nonVax accounts for the (inefficient) administration of vaccines to exposed, infectious and hospitalised individuals
%nonVax approximates S+E+I+H+R (unvaccinated) and D (partially)
%S (or R) ./nonVax is approximately 1 (or 0) when prevalence is low but is closer to 0 (or 1) when prevalence is high
%nonVax is non-zero as long as uptake is less than 100%

if t>=vx.end
    v1rates=zeros(ntot,1);
    v1rater=zeros(ntot,1);
    Vdot=   zeros(ntot,1);
    
elseif t>=vx.startp5
    v1rates=vx.ratep5period.*S./nonVax;
    v1rater=vx.ratep5period.*R./nonVax;
    Vdot=   vx.ratep5period;
       
elseif t>=vx.startp4
    v1rates=vx.ratep4period.*S./nonVax;
    v1rater=vx.ratep4period.*R./nonVax;
    Vdot=   vx.ratep4period;
    
elseif t>=vx.startp3
    v1rates=vx.ratep3period.*S./nonVax;
    v1rater=vx.ratep3period.*R./nonVax;
    Vdot=   vx.ratep3period;
    
elseif t>=vx.startp2
    v1rates=vx.ratep2period.*S./nonVax;
    v1rater=vx.ratep2period.*R./nonVax;
    Vdot=   vx.ratep2period;
    
elseif t>=vx.startp1
    v1rates=vx.ratep1period.*S./nonVax;
    v1rater=vx.ratep1period.*R./nonVax;
    Vdot=   vx.ratep1period;
    
else
    v1rates=zeros(ntot,1);
    v1rater=zeros(ntot,1);
    Vdot=   zeros(ntot,1);
    
end
%
if t<472%T_delta_introduced 16th April
    vxscv1=vx.scv1;
elseif t>548%t_delta_only 1st July
    vxscv1=vx.scv1b;
else
    vxscv1=vx.scv1b+(t-472)/76*(vx.scv1-vx.scv1b);
end

%% EQUATIONS:

Sdot=       -S.*(foi+seed)*redSus  +pr.nu.*R   -v1rates                +vx.nuv1.*Sv1;
Edot=        S.*(foi+seed)*redSus  +Shv1.*foi  -pr.sigma.*E;

Inadot=     (1-p3)      *(1-pr.p1)              *pr.sigma   .*E         -pr.g1*(1+pr.odds)          .*Ina;
Isadot=     p3          *(1-pr.p1)              *pr.sigma   .*E         -pr.g1*(1+pr.odds)          .*Isa;
Insdot=     (1-p4)      *pr.p1                  *pr.sigma   .*E         -(pr.g2+pr.h)               .*Ins;
Issdot=     p4          *pr.p1                  *pr.sigma   .*E         -(pr.g2+pr.h+pr.q2)         .*Iss;

Hdot=       pr.h.*(Ins+Iss)                                             -(pr.g3+pr.mu).*(min(occ,pr.Hmax)/occ).*H           -(pr.g3_oc+pr.mu_oc).*(max(0,occ-pr.Hmax)/occ).*H;
Rdot=       pr.g1.*(1+pr.odds).*(Ina+Isa)+pr.g2.*(Ins+Iss)              + pr.g3.*(min(occ,pr.Hmax)/occ).*H                  + pr.g3_oc.*(max(0,occ-pr.Hmax)/occ).*H         -pr.nu.*R   -v1rater;

Shv1dot=    v1rates     -vx.hrv1*Shv1   -Shv1.*foi;
Sv1dot=                  vx.hrv1*Shv1   -Sv1.*(1-vxscv1).*foi  -vx.nuv1.*Sv1;  %+pr.nu.*Rv1;   
Ev1dot=                                  Sv1.*(1-vxscv1).*foi  -pr.sigma.*Ev1;

Inav1dot=   (1-p3)      *(1-pr.p1*(1-vx.p1v1))  *pr.sigma   .*Ev1       -pr.g1*(1+pr.odds)          .*Inav1;
Isav1dot=   p3          *(1-pr.p1*(1-vx.p1v1))  *pr.sigma   .*Ev1       -pr.g1*(1+pr.odds)          .*Isav1;
Insv1dot=   (1-p4)      *(1-vx.p1v1)*pr.p1      *pr.sigma   .*Ev1       -(vx.g2_v1+vx.h_v1)         .*Insv1;
Issv1dot=   p4          *(1-vx.p1v1)*pr.p1      *pr.sigma   .*Ev1       -(vx.g2_v1+vx.h_v1+pr.q2)   .*Issv1;

Hv1dot=     vx.h_v1.*(Insv1+Issv1)                                      -(vx.g3_v1+vx.mu_v1).*(min(occ,pr.Hmax)/occ).*Hv1   -(vx.g3_ocv1+vx.mu_ocv1).*(max(0,occ-pr.Hmax)/occ).*Hv1;
Rv1dot=     pr.g1.*(1+pr.odds).*(Inav1+Isav1)   +vx.g2_v1.*(Insv1+Issv1)+ vx.g3_v1.*(min(occ,pr.Hmax)/occ).*Hv1             + vx.g3_ocv1.*(max(0,occ-pr.Hmax)/occ).*Hv1     +v1rater;%-pr.nu.*Rv1;

DEdot=      pr.mu.*(min(occ,pr.Hmax)/occ).*H    +pr.mu_oc .*(max(0,occ-pr.Hmax)/occ).*H     +vx.mu_v1.*(min(occ,pr.Hmax)/occ).*Hv1  +vx.mu_ocv1.*(max(0,occ-pr.Hmax)/occ).*Hv1;     

incidence=S.*(foi+seed)*redSus+Shv1.*foi+Sv1.*(1-vxscv1).*foi;

%% OUTPUT:

f= [Sdot;Edot;...
    Inadot;Isadot;Insdot;Issdot;...
    Hdot;Rdot;...
    Shv1dot;Sv1dot;Ev1dot;...
    Inav1dot;Isav1dot;Insv1dot;Issv1dot;...
    Hv1dot;Rv1dot;...
    DEdot;Vdot;incidence];

g=sum(pr.h.*(Ins+Iss)  + vx.h_v1.*(Insv1+Issv1)  );%Hin
%g=S.*(foi+seed)*redSus+Shv1.*foi+Sv1.*(1-vxscv1).*foi; %Incidence

end

%%

%Stochastic variant - needs update to C19 flowchart
function f=stochSim(y,beta,gamma,n,ntot,NN,N0,D,seed,phi1,phi2,tau,alpha)
%Still flu-like/SIR****
%Feed in mu if required
factor=6;
tend=360*factor; beta=beta/factor; gamma=gamma/factor;
Vec=zeros(ntot,tend);
%
S=y(1:ntot);
I=y(ntot+1:2*ntot);
R=y(2*ntot+1:end);
i=1;
threshold=30;%Number of time-steps for necessary simulation/seed
while i<tend && (i<threshold || sum(I)>0)%At least 30 time-steps
phi=1;%phi1-phi2*cos(pi*i*f1/180);%Just seed here
Sout=1-exp(-phi*(beta*(D*(I./N0).^alpha)+seed*heaviside(threshold-i)));%+mu*R;.^alpha
Sout(Sout>1)=1;
Sout=binornd(S,Sout); Sout(S==0)=0;
S=S-Sout; I=I+Sout;
%
Iout=1-exp(-gamma);
Iout=binornd(I,Iout);
I=I-Iout; R=R+Iout;
Vec(:,i)=I;
i=i+1;
end
f=R;
if sum(isnan(R))>0
    fuck=1;
    print('Somenting is NaN')
end
end

%%

function [value,isterminal,direction]=changeinbehave(t,y,ntot,pr,vx)

    H=      y(6*ntot+1:7*ntot);
    Hv1=    y(15*ntot+1:16*ntot);
    occ=    max(1,sum(H+Hv1));

    value      = min(occ-pr.Hmax,0)+min(t-pr.Tm-0.1,0);
    direction  = 0;
    isterminal = ones(size(value));
    
end

function [value,isterminal,direction]=lockdown(t,y,ntot,pr,vx)

    Ins=    y(4*ntot+1:5*ntot);
    Iss=    y(5*ntot+1:6*ntot);
    H=      y(6*ntot+1:7*ntot);
    Insv1=  y(13*ntot+1:14*ntot);
    Issv1=  y(14*ntot+1:15*ntot);
    Hv1=    y(15*ntot+1:16*ntot);
    occ=    max(1,sum(H+Hv1)); 
    Hdot=   pr.h.*(Ins+Iss)         -(pr.g3+pr.mu).*(min(occ,pr.Hmax)/occ).*H           -(pr.g3_oc+pr.mu_oc).*(max(0,occ-pr.Hmax)/occ).*H;
    Hv1dot= vx.h_v1.*(Insv1+Issv1)  -(vx.g3_v1+vx.mu_v1).*(min(occ,pr.Hmax)/occ).*Hv1   -(vx.g3_ocv1+vx.mu_ocv1).*(max(0,occ-pr.Hmax)/occ).*Hv1;
    occdot= sum(Hdot+Hv1dot);
     
    r=occdot/occ;
    Tcap=t+log(pr.Hmax/occ)/r;
    Tcap=Tcap-7;

    value      = [min(t-Tcap,0)+min(t-pr.Tm-0.1,0),min(occ-0.95*pr.Hmax,0)+min(t-pr.Tm-0.1,0)];
    direction  = [1,1];
    if r>0.025
        isterminal(1) = 1;
    else
        isterminal(1) = 0;
    end
    isterminal(2) = 1;
    
end

function [value,isterminal,direction]=reopen(t,y,ntot,pr,vx)

    H=      y(6*ntot+1:7*ntot);
    Hv1=    y(15*ntot+1:16*ntot);
    occ=    max(1,sum(H+Hv1));

    value      = occ-pr.thl;
    direction  = -1;
    isterminal = ones(size(value));
    
end

%%

function f=plotEpi(tout,Y,H,n,dodiff,tvec)
solvetype=2;
tend=tvec(end);%720;%For plot only;
cmap=lines(7);
if solvetype==2
    figure
    fs=10; lw=2;
    Hall=sum(H,2);
    maxY=max(max(Hall));
    maxY=1e3*ceil(maxY/1e3);
    hold on
    for i=2:length(tvec)-1
        plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
    end
    plot(tout,Hall,'-','linewidth',lw,'color',cmap(2,:));%'--', (i,:)
    xlabel('Time','FontSize',fs);
    ylabel('Population','FontSize',fs);
    set(gca,'FontSize',fs);
    axis([tvec(1),tend,0,maxY])
    grid on
    grid minor
    box on
    hold off
end
end

function f=yofh(H)
%Doesn't allow for negative H
hmin=5e3;
hmax=5e4;
ymin=.3;
if H<hmin
    f=1;
elseif H<hmax
    f=1-(H-hmin)/(hmax-hmin)*(1-ymin);
else
    f=ymin;
end
end

function f=rho2ofh(H,pr)
%%**
f=pr.m1+(1-pr.m1)*pr.L1/(1+exp(dot(pr.k1,H)-pr.H01'));
%f=pr.L1;

%f=pr.m1+(1-pr.m1)*pr.L1/(1+exp(dot(pr.k1,H)-pr.H01));
%)
%f=min(1,1-pr.L1*pr.k1*(H-pr.H01)/(1+pr.L1*(H-pr.H01)));
%f=pr.L1*exp(-pr.k1*H);
%f=pr.L1+pr.k1*H;
end

function f=beFeedback(H,pr)
f=pr.m1+(1-pr.m1)*pr.L1/(1+exp(dot(pr.k1,H)-pr.H01'));
end