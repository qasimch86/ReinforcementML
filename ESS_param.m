function par=ESS_param()
global k1 k2 k3 k4 re1 re2 Rm Rd Pdiv Pdeath K div
k1=1.5*10^4;%9.16;%10^7;%
k2=log(2);
k3=0.75;%k3=[1.4 1.9 2.1]
k4=log(2);
% re1=0;%0.125 
% re2=0;%0.125 
% Rm=0.5;%0.75
% Rd=1-Rm;
Pdiv=1500;%+div*15;
Pdeath=600;%+div*5;
K=2500;%+div*20;