function[IRCF,LCF,MLCF,CCF,MCI,CTT,CA,ST,VCF] =PF()
a=75;
b=100;
d=12;
e=270;
f=340;
g=50;
h=100;
i=303;
j=333;

IRCF=a + (b-a)*rand();
LCF=a + (b-a)*rand();
MLCF=a + (b-a)*rand();
CCF=a + (b-a)*rand();
MCI=randi(d);
CTT= e + (f-e)*rand();
CA = g + (h-g)*rand();
ST= i + (j-i)*rand();
VCF = a + (b-a)*rand();
end