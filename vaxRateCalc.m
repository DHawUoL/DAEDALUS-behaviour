%{
4.7e7 1st doses by 31st July
90% coverage in old people

95% of 11814000 between 10th Dec 14th Feb
11223300 between days 345 and 411
=170050/day

38378200 total working age
4.7e7-11223300=35776700

35776700/(579-411)=2.1296e+05

2.1296e+05*30550016/38378200=1.6952e+05 working rate

2.1296e+05*7828184/38378200=4.3438e+04 non-working rate

So:

startp1: 345
startp2: 411
startp3: 579
startp4: 580
startp5: 581

aratep1=[0;0;0;170050]
aratep2=[0;0;2.1296e+05;0];
aratep3=zeros(4,1);
aratep4=zeros(4,1);
aratep5=zeros(4,1);

%}