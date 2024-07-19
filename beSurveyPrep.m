%A=beTable(2+[5,8,9,10],2:end);
A=beTable(2+[5,8,10],2:end);
[N,M]=size(A);
out=zeros(N,M);
for i = 1 : N
   for j = 1 : M
        val = A{i,j};
        out(i,j)=str2double(string(val));
    end
end

crudeAv=[2,7,12,16,21,24,28,32,36,40,44,48,53,57,61,66,70];
lc=length(crudeAv);
avInds=zeros(size(out,1),lc-1);
for i=1:lc-1
    diff=crudeAv(i+1)-crudeAv(i);
    avInds(:,i)=sum(out(:,crudeAv(i):(crudeAv(i+1)-1)),2)/diff;
end
