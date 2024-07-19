%{
If PCR:
coeffs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
If hospital coccupancy etc:

%%BH flag
Respect length(coeff) as number of input dimensions to feedback
L1, K1, H01
beProps - find

Xit - single sector fully open
OR input to logistic function
%}

[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(data,R0,del,forFeedback,coeff,beProps,alphaIn);
%Second argument, "numInt", gone
[pr,be,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=bePrepCovid19(dataUK1,3,ones(1,size(X,2)-2),[.8,.01,.01,.01,0],ones(1,3)',zeros(5,1),[.3,.3,.3]);
%Function arguments chamged - Sim also
[f,g,tvec,rhohat]=beRunCovid19(pr,be,vx,n,NN,NNbar,beta,ones(1,length(tvec)-1),tvec,1,dataUK1);

%In sim code:
if i==1
                propsBi=zeros(5,1);
            else
                %{
                Hx1=sum(Hout(end,:));
                Hx2=sum(HnewAll(end,:));
                Hx3=Hx1-sum(Hout(end-13,:));
                propsBi=beFeedback2([Hx1,Hx2,Hx3],pr)*ones(5,1);%ones(5,1);%
                %}
                propsBi=be.BiFirstFit*ones(5,1);
end

[poptim,Ypred,delta,resnorm]=beFitEpiStart(dataAdmissions,ones(1,size(X,2)-2),dataUK1,[2.8,-70,85,0,1],ones(1,size(X,2)-2),ones(1,3)');
1.9254 -195.9958   93.5537    0.2050    0.5319

[poptim,Ypred,delta,resnorm]=beFitEpi(dataAdmissions,ones(1,size(X,2)-2),dataUK1,[],ones(1,size(X,2)-2),ones(1,3)');

%X for makeDs, Xfull for feedback
x=stringToJuly21;
x=x([1,2,4:19],1)';
[poptim,Ypred,delta,resnorm]=beFitEpi(dataAdmissions,ones(1,size(X,2)-2),dataUK1,[2.8,-70,85,0,1],x,ones(1,2)');
0.3810    0.6130   -0.0256   -4.7365  -13.3980

x=stringToJuly21;
x=x([1,3:19],[2,5,11])';
[poptim,Ypred,delta,resnorm]=beFitEpi(dataAdmissions,ones(1,size(X,2)-2),dataUK1,[0.5    1   0  -4   -6  -50],x,ones(1,3)');
%Indices 1&11, 2 alphas:
0.2638    0.8000    0.5893   -7.5405  -10.0000  -50.0000
0.5174    0.8623    0.8191   -8.2332   -5.2586  -77.5589
%Indices 5&11, 2 alphas:
0.5174    0.8623    0.8191   -8.2332   -5.2586  -77.5589
%Indices 2,5,11, 2 alphas:
[poptim,Ypred,delta,resnorm]=beFitEpi(dataAdmissions,ones(1,size(X,2)-2),dataUK1,[0.5174,0.8623,0.8191,2,-8.2332,-5.2586,77.5589],x,ones(1,3)');
%Indices 2,11, 2 alphas:
x=stringToJuly21;
x=x([1,3:19],[2,11])';
[poptim,Ypred,delta,resnorm]=beFitEpi(dataAdmissions,ones(1,size(X,2)-2),dataUK1,[0.5168    0.8624    0.8196  -20.0000   -5.2586  -75.0000],x,ones(1,2)');

%fitEpiFull
%if i>4 add modifier ********, i<3 line 191, alpha 2 from i=9, 
%i>4 in Bayesian plot
x=stringToJuly21;
x=x(:,[2,11])';
[poptim,Ypred,delta,resnorm]=beFitEpiFull(dataAdmissions,ones(1,size(X,2)-2),dataUK1,[2,-50,0.5,0.8,0.8,-20,-5.2,-75],x,ones(1,2)');
x=stringToJuly21;
x=x(:,[1,2,5])';%R0 fixed/range; 1st 3 shouldnt matter
[poptim,Ypred,delta,resnorm]=beFitEpiFull(dataAdmissions,ones(1,size(X,2)-2),dataUK1,[2.8,-68,0.3,0.8,0.8,-1,-3,-2,-50],x,ones(1,3)');
%Fix R0=2.8, t0=-68:
[poptim,Ypred,delta,resnorm]=beFitEpiFull(dataAdmissions,ones(1,size(X,2)-2),dataUK1,[0.3,0.8,0.8,-1,-3,-2,-50],x,ones(1,3)');
0.2762    0.6864    0.9279   -0.4202   -1.2592   -2.0000 -113.1649
0.2507    0.6148    0.7750   -9.1953   -9.4777   -9.9271 -100.7279
%Fix R0=2.75, t0=-70:
0.3065    0.75    0.7154   -3   -4  -14  -80
0.3065    0.65    0.7154   -3   -4  -14  -80 - manual
0.0023    0.5672    0.5679  -19.0766   -4.6897   -6.4059  -49.7466 %String, gov resp, school closures

%Bayesian:
[xsto, outsto, history, accept_rate,covmat]=beFitEpiBayesian(dataAdmissions,ones(1,size(X,2)-2),dataUK1,poptim,x,ones(1,3)');
bePlotBayesianFitEpi(dataAdmissions,ones(1,size(X,2)-2),dataAdmissions,ones(1,size(X,2)-2),xsto,x,ones(1,3)',poptim);

x=stringToJuly21';
bePlotBayesianFitEpi(dataAdmissions,ones(1,size(x,2)),dataUK1,xsto1b,x,ones(1,1)',poptim1);
x=stringToJan22(:,1)';
bePlotBayesianFitEpi(dataAdmissionsExtended,ones(1,size(x,2)),dataUK1,xsto1b,x,ones(1,1)',poptim1);

%Fitting individual p's:
[poptim,Ypred,delta,resnorm]=beFitEpiFull(dataAdmissionsExtended,ones(1,size(x,2)),dataUK1,[.2,.5*ones(1,32)],x,1);
%One alpha:
[0.2645    0.6990    0.6577    0.6160    0.6519    0.4436    0.25    0.25    0.6    0.25    0.4810    0.7604    0.4027    0.3970    0.2352    0.2   0.45    0.55    0.35    0.4    0.45    0.1    0.3808]
%Upper bounds 7,8,10=.35:
[0.2647    0.6998    0.6565    0.6149    0.6514    0.4443    0.2513    0.2503    0.6019    0.2493    0.4799    0.7611    0.1623    0.4360    0.3001    0.2526    0.5006 0.5089    0.4928    0.3624    0.5000    0.1582    0.3256];
%{
Notes for talk:
Take out March, switch at date calibrated
%Fix R0, switch at end of March
Epi info not enough - c.f. November lockdown 
%}