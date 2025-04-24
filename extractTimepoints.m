function xout=extractTimepoints(x1,x2,x3)

t1=x1(:,1); x1=x1(:,2);
t2=x2(:,1); x2=x2(:,2);

t3=x3(:,1); x3=x3(:,2);

points=unique([t1;t2;t3]);
lp=length(points);
xout=[sort(points),zeros(lp,3)];%t, x1, x2, x3

x1hold=x1(1);
x2hold=x2(1);

x3hold=x3(1);

for i=1:lp
    pointi=points(i);
    f1=find(t1==pointi);
    f2=find(t2==pointi);

    f3=find(t3==pointi);

    if isempty(f1)==0
        x1i=x1(f1);
        xout(i,2)=x1i;
        x1hold=x1i;
    else
        xout(i,2)=x1hold;
    end
    if isempty(f2)==0
        x2i=x2(f2);
        xout(i,3)=x2i;
        x2hold=x2i;
    else
        xout(i,3)=x2hold;
    end

    if isempty(f3)==0
        x3i=x3(f3);
        xout(i,4)=x3i;
        x3hold=x3i;
    else
        xout(i,4)=x3hold;
    end

end
