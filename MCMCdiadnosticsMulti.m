function f=MCMCdiadnosticsMulti(cell)
nChains=length(cell);
parameter=1;

accept=zeros(1,nChains);
newIC=zeros(nChains,3);
figure
hold on
for i=1:nChains
    xstoi=cell{i}.xsto;
    plot(xstoi(:,parameter))
    newIC(i,:)=xstoi(end,:);
    %accept(i)=cell{i}.accept_rate;
end
%disp(accept)

f=newIC;

