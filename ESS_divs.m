function Fdiv=ESS_divs(time,rm)
% global re1 re2 Rm Rd Pdiv npas P Pd div_t Pdeath
global Pdiv npas Rm Rd re1 re2 Pdeath Altrusm div_t remax P
n=1;
P=zeros(2,npas);
Pd=zeros(2,npas);
%% Finding D0
P(1,1)=Pdiv*Rd;
P(2,1)=0;
div=0;
parent=0;
daugh=0;val=0;
ESS_param();
div=0;RE=0;re1=1;re2=re1;
flag=1;
%% RK Method
        sp=ESS_rhs(P(:,1));
        In_sp(1)=sp(1);Dm_sp(1)=sp(2);
%         RE(1,1)=re1;
for k=1:npas-1
    %% RK Method
    if P(1,k)<Pdiv
        P(:,k+1)=ESS_desol(P(:,k));
        sp=ESS_rhs(P(:,k+1));
        In_sp(k+1)=sp(1);Dm_sp(k+1)=sp(2);
    else
        div=div+1;
        div_t(div+1)=k;
%         if flag==1
            re1=1-((div)^2/((div)^2+1));%(P(2,div_t(div)+1)/P(2,div_t(div+1))));%
            re1=re1+2*(Altrusm-1/2)*((1-Altrusm)*re1+Altrusm*(remax-re1));
%         end
%         DIV_P=protein_div(P(:,k),0);
%         Q=ESSdivs_slv(DIV_P(1:2)');
%         while Q(1)<0%Pdiv
%             re1=re1-0.1*re1;re2=re1;
%             if re1<=0.01
%                 re1=0;
%                 break;
%             end
%             flag=0;
%             DIV_P=protein_div(P(:,k),0);
%             Q=ESSdivs_slv(DIV_P(1:2)');
%         end
        if re1>1
            re1=1;
        end
        if re1<0.01
            re1=0;
        end
        re2=re1;
        RE(div,1)=re1;
        %% for mother
        DIV_P=protein_div(P(:,k),1);
        P(1,k+1)=DIV_P(1);
        P(2,k+1)=DIV_P(2);
        sp=ESS_rhs(P(:,k+1));
        In_sp(k+1)=sp(1);Dm_sp(k+1)=sp(2);
        %% for daughter
        Pd(1,k+1)=DIV_P(3);
        Pd(2,k+1)=DIV_P(4);
        parent(div)=P(2,k+1)/(P(2,k)); %(P(2,k+1))/Psum;%P(2,k)/P(2,k+1);%P(2,k+1)/P(1,k+1);%P(2,k+1)/sum;%((P(2,k)*Rm)/P(2,k+1))/((P(1,k)*Rm)/P(1,k+1));
        daugh(div)=Rm/Rd*Pd(2,k+1)/(P(2,k)); %Pd(2,k+1)/Pdsum;%Pd(2,k+1)/Pd(1,k+1);%Pd(2,k+1)/Sum;%Pd(2,k+1)/Pdsum;%((Pd(2,k)*Rm)/Pd(2,k+1))/((Pd(1,k)*Rd)/Pd(1,k+1));
    end
end
auto_P(1:div,1)=parent;
auto_D(1:div,1)=daugh;
Fdiv=[auto_P;auto_D;RE;div];