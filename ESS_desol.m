function P=ESS_desol(P0)
global delt
K1 = ESS_rhs(P0);
K2 = ESS_rhs(P0 + delt*K1./2);
K3 = ESS_rhs(P0 + delt*K2./2);
K4 = ESS_rhs(P0 + delt*K3);
P=P0+ delt*(K1 + 2*K2 + 2*K3 + K4)/6;