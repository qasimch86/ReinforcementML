function DIV=protein_div(P0,flag)
global Rm Rd re1 re2
        P1(1)=P0(1)*Rm-P0(2)*Rd*re1;
        P1(2)=P0(2)*Rm+P0(2)*Rd*re2;
        if flag==1
            Pd(1)=P0(1)*Rd+P0(2)*Rd*re1;
            Pd(2)=P0(2)*Rd-P0(2)*Rd*re2;
        else
            Pd(1)=0;
            Pd(2)=0;
        end
        DIV=[P1 Pd];