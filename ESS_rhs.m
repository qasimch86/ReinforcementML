%  user defined function for Standard RK Method
function f=ESS_rhs(P)
global k1 k2 k3 k4 re Rm Rd K count Pdiv

I0=P(1);
D0=P(2);

f=[k1*(1-(I0+D0)/K)-k2*I0-k3*I0;%Pintact
   k3*I0-k4*D0];
% f=[k1/(100+I0+D0)-k2*I0-k3*I0;%Pintact
%    k3*I0-k4*D0];
   %r1*Sc/Vc*c(1);%Pdamage