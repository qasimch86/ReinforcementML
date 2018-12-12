%% Initial values for Pdam and Pint
% clear all
function ESS_genval()
tic
global k1 k2 k3 k4 re1 re2 Rm Rd Pdiv Pdeath K count delt npas tmax
% tmax=1;
% npas=301;
% delt=(tmax/(npas-1));
time=0:delt:tmax;
% Rm=0.79;Rd=1-Rm;
ESS_param();
IC=[Pdiv*Rd 0];
n=1/4; nre=0;
%% Initial values Ranges
MtIn_sz=1;%round(Pdeath*Rd/5);
DgIn_sz=1;%round(Pdeath*Rd/5);
MtDm_sz=1;%round(Pdeath*Rm/5);
DgDm_sz=1;%round(Pdeath*Rd/5);
min_stsz=DgDm_sz;
MtIn_givl=round(Pdiv*Rm-Pdeath*Rd):MtIn_sz:round(Pdiv*Rm)+10;%:1000:1185+15;
MtIn_givl_len=length(MtIn_givl);
DgIn_givl=round(Pdiv*Rd):DgIn_sz:round(Pdiv*Rd+Pdeath*Rd);%315:420;% Daugther intact initial value
DgIn_givl_len=length(DgIn_givl);
MtDm_givl=0:MtDm_sz:Pdeath*Rm+Pdeath*Rd;%0:500;% Mother damage initial value
MtDm_givl_len=length(MtDm_givl);
DgDm_givl=0:DgDm_sz:round(Pdeath*Rd);%0:105;
DgDm_givl_len=length(DgDm_givl);
%% Initialization
Mt_div(1:MtIn_givl_len,1:MtDm_givl_len)=0;% Mother time to division
Dg_div(1:DgIn_givl_len,1:DgDm_givl_len)=0;% Daughter time to division
MtIn_gfvl(1:MtIn_givl_len,1:MtDm_givl_len)=0;% Mother Intact final value
DgIn_gfvl(1:DgIn_givl_len,1:DgDm_givl_len)=0;% Daughter intact final value
MtDm_gfvl(1:MtIn_givl_len,1:MtDm_givl_len)=0;% Mother damage final value
DgDm_gfvl(1:DgIn_givl_len,1:DgDm_givl_len)=0;% Daughter damage final value
%% For Initial and final values for pdam from 0 to Pdeath
iIv=0;
for Ival=MtIn_givl%1000:1125
    iIv=iIv+1;iDv=0;
    for Dval=MtDm_givl%0:375
        iDv=iDv+1;
        M=zeros(2,npas);
        M(:,1)=[Ival;Dval];
        for k=1:npas-1
        %% RK Method
            if M(1,k)<Pdiv
                M(:,k+1)=ESS_desol(M(:,k));
            else
                MtIn_gfvl(iIv,iDv)=M(1,k);
                MtDm_gfvl(iIv,iDv)=M(2,k);
                Mt_div(iIv,iDv)=k;
                break;
            end
            if M(2,k)>Pdeath
                break;
            end
        end 
    end
end
iIv=0;
for Ival=DgIn_givl%375:500
    iIv=iIv+1;iDv=0;
    for Dval=DgDm_givl%0:125
        iDv=iDv+1;
        D=zeros(2,npas);
        D(:,1)=[Ival;Dval];
        for k=1:npas-1
        %% RK Method
            if D(1,k)<Pdiv
               D(:,k+1)=ESS_desol(D(:,k));
            else
               DgIn_gfvl(iIv,iDv)=D(1,k);%Daughter intact final value
               DgDm_gfvl(iIv,iDv)=D(2,k);%Daughter damage final value
               Dg_div(iIv,iDv)=k;%Daughter division time duration
               break;
            end
            if D(2,k)>Pdeath
               break;
            end
        end
    end
end
save('ESS_genval.mat','MtIn_givl','MtIn_gfvl','MtDm_givl','MtDm_gfvl','DgIn_givl','DgIn_gfvl','DgDm_givl','DgDm_gfvl','Mt_div','Dg_div');
end