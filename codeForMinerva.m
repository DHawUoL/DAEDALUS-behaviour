yinfinal3=Keydatestablev2;
tvec3=yinfinal3(:,1)';
yinfinal3=yinfinal3(:,2:3)';
yinfinal3=[zeros(2,3),yinfinal3(:,3:end)];
yinfinal3(2,:)=yinfinal3(2,:)/100;


tvec3(1)=1;
tvec3(end+1)=1091;
tdiff=diff(tvec3);
toPlot=[repelem(yinfinal3(2,2:end)',tdiff')];

finalData3=finalData2;
finalData3{:,4}=-toPlot;