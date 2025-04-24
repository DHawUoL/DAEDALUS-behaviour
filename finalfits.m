
[ystofinal2e,~,accept_rate]=beFitEpiBayesian(dataAdmissions,ones(1,size(yinfinal2,2)),dataUK2,[ystofinal2d(end,:)],yinfinal2,ones(1,3)');

[ysto,~,accept_rate]=beFitEpiBayesian_alpham(dataAdmissions,ones(1,size(yinfinal2,2)),dataUK2,[ystofinal2e(end,:)],yinfinal2,ones(1,3)');
a=.2293;%Also hard-coded in sim2fit
b=.6942;
ystom=[ysto(:,1),a*ysto(:,1)+b,ysto(:,2:end)];

[ysto2,~,accept_rate]=beFitEpiBayesian_alpham(dataAdmissions,ones(1,size(yinfinal2,2)),dataUK2,[ysto(end,:)],yinfinal2,ones(1,3)');
%-different x0

[ysto2,~,~,accept_rate,covmat]=beFitEpiBayesian_alpham(dataAdmissions,ones(1,size(yinfinal2,2)),dataUK2,[ysto2(end,:)],yinfinal2,ones(1,3)');
%-sigma=0.5
%-this is c_plot - use this one

[ysto3,~,~,accept_rate,covmat]=beFitEpiBayesian_alpham(dataAdmissions,ones(1,size(yinfinal2,2)),dataUK2,[ysto2(end,:)],yinfinal2,ones(1,3)');
%-sigma=0.2

%xsto1block
%xsto2block - differnt sigmas, started from end of previous - identiftying linear relationships