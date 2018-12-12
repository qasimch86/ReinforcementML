%% This program finds population on dynamic retention values
clear all
tic
global Rm Rd Pdeath delt npas tmax Pdiv re1 re2 Altrusm
tmax=1;
dtmax=300*tmax;
scl_fct=1;
npas=round(scl_fct*dtmax)+1;
delt=(tmax/(npas-1));
time=0:delt:tmax;
Rm=0.79;Rd=1-Rm;
re2=1;re1=re2;
rm=[Rm];len_rm=length(rm);
%%%%%%%%%%%%%%%%%%% General Values
if exist('ESS_genval.mat')
    load ESS_genval.mat;
else
    ESS_genval();
    load ESS_genval.mat
end
%%%%%%%%%%%%%%%%%%%
u=1;
Tot_Dm=0;
TotSizewise_payoff(:,1)=[10 20];
Altrusm=[1;1];
while abs(diff(TotSizewise_payoff(:,u)))>.1
ESS_param();
% while abs((diff(Tot_Dm))>0.1
u=u+1;
for st=1:2
re_vec(u,st,1:20)=0;
    if st==1
        Strgy1=ESS_dist(time,rm);%0.875:-0.2:0;%[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];%[0 0.3 0.6 1];%
        div1=Strgy1(end);
        re_vec(u,st,1:div1)=Strgy1(2*div1+1:3*div1);
%         re_vec=Altrusm(1);
    elseif st==2
        Strgy2=ESS_divs(time,rm);
        div2=Strgy2(end);
        re_vec(u,st,1:div2)=Strgy2(2*div2+1:3*div2);
    end
IC=[Pdiv*Rd 0];
re_len=length(re_vec(u,st,:));
%% Retention dependence upon anticipated number of possible divisions
    N_Dg(1)=1;N_Mt=0;
    Dg_flg=1;Mt_flg=0;
    nm_Dg_Mt=0;div_sm=0;
    MtDiv_nm=zeros(npas+100,npas);
    Dg_tdiv=0;DgIn_fvl=0;DgDm_fvl=0;
    Mt_tdiv=0;MtIn_fvl=0;MtDm_fvl=0;
    MtDiv_nmOld=0;Mt_ivl=0;Dg_ivl=0;
    Dg_ivl(1,1:2)=IC;nm_Dgdiv=0;
for k=1:npas-1
    Time = k
    Mt_ivl_N=0;
    Dg_ivl_N=0;
    %% for daughter
  if Dg_flg==1
    for Dnum=1:N_Dg
        Dg_ivl_N(1:2)=round(Dg_ivl(Dnum,1:2));%Intact and damge comp. of daughter cell 
        X=find((DgIn_givl==Dg_ivl_N(1)));
        Y=find((DgDm_givl==Dg_ivl_N(2)));
      if isempty(X)
          [x0,x1]=min(abs(DgIn_givl-Dg_ivl_N(1)));
          X=find((DgIn_givl==DgIn_givl(x1)));%find intact initial value from general daughter values
          sprintf('Daughter Values changed: from Intact = %.0f to Intact = %.0f',Dg_ivl_N(1),X)
      end
      if isempty(Y)
          [y0,y1]=min(abs(DgDm_givl-Dg_ivl_N(2)));
          Y=find((DgDm_givl==DgDm_givl(y1)));%find intact initial value from general daughter values
          sprintf('Daughter Values changed: from Damage = %.0f to Damage = %.0f',Dg_ivl_N(2),Y)
      end
        Dg_tdiv(k,Dnum)=round(scl_fct*Dg_div(X,Y))+k-1;% Daughter division time
        Dg_tdiv_sm=sum(sum(ismember(Dg_tdiv,Dg_tdiv(k,Dnum))));
        DgIn_fvl(Dg_tdiv(k,Dnum),Dg_tdiv_sm)=DgIn_gfvl(X,Y);% Daughter intact final value
        DgDm_fvl(Dg_tdiv(k,Dnum),Dg_tdiv_sm)=DgDm_gfvl(X,Y);% Daughter damage final value
    end
%     DgIn_fvl_vec(k,1:max(Dg_tdiv(k,:)),1:Dg_tdiv_sm)=DgIn_fvl(1:max(Dg_tdiv(k,:)),1:Dg_tdiv_sm);
%     DgDm_fvl_vec(k,1:max(Dg_tdiv(k,:)),1:Dg_tdiv_sm)=DgDm_fvl(1:max(Dg_tdiv(k,:)),1:Dg_tdiv_sm);
    Dg_flg=0;N_Dg=0;Dg_ivl=0;
  end
    %% for mother
  if Mt_flg==1
    for Mnum=1:N_Mt
        Mt_ivl_N(1:2)=round(Mt_ivl(Mnum,1:2));%Intact and damge comp. of daughter cell 
        X=find((MtIn_givl==Mt_ivl_N(1)));
        Y=find((MtDm_givl==Mt_ivl_N(2)));
      if isempty(X)
          [x0,x1]=min(abs(MtIn_givl-Mt_ivl_N(1)));
          X=find((MtIn_givl==MtIn_givl(x1)));
          sprintf('Mother Values changed: from Intact = %.0f to Intact = %.0f',Mt_ivl_N(1),MtIn_givl(x1))
      end
      if isempty(Y)
          [y0,y1]=min(abs(MtDm_givl-Mt_ivl_N(2)));
          Y=find((MtDm_givl==MtDm_givl(y1)));
          sprintf('Mother Values changed: from Damage = %.0f to Damage = %.0f',Mt_ivl_N(2),MtDm_givl(y1))
      end
        Mt_ivl_srch(1)=X;%find intact initial value from general Mother values
        Mt_ivl_srch(2)=Y;%find intact initial value from general Mother values
        Mt_tdiv(k,Mnum)=round(scl_fct*Mt_div(Mt_ivl_srch(1),Mt_ivl_srch(2)))+k-1;% Mother division time
        Mt_tdiv_sm=sum(sum(ismember(Mt_tdiv,Mt_tdiv(k,Mnum))));
      if Mt_tdiv(k,Mnum)>k
        MtIn_fvl(Mt_tdiv(k,Mnum),Mt_tdiv_sm)=MtIn_gfvl(Mt_ivl_srch(1),Mt_ivl_srch(2));% Mother intact final value
        MtDm_fvl(Mt_tdiv(k,Mnum),Mt_tdiv_sm)=MtDm_gfvl(Mt_ivl_srch(1),Mt_ivl_srch(2));% Mother damage final value
        if Mnum<=nm_Dgdiv
            MtDiv_nm(Mt_tdiv(k,Mnum),Mt_tdiv_sm)=1;
        elseif Mnum>nm_Dgdiv
            MtDiv_nm(Mt_tdiv(k,Mnum),Mt_tdiv_sm)=MtDiv_nmOld(Mnum)+1;
        end
%         if Mt_tdiv(k,Mnum)>=npas-1
%             MtAlive_ivl(n,:)=[Mt_tdiv(k,Mnum); Mt_ivl_srch];
%         end
      else
        MtIn_fvl(Mt_tdiv(k,Mnum),Mt_tdiv_sm)=0;% probably when Mt_div = 0
        MtDm_fvl(Mt_tdiv(k,Mnum),Mt_tdiv_sm)=0;% probably when Mt_div = 0
        MtDiv_nm(Mt_tdiv(k,Mnum),Mt_tdiv_sm)=0;
        Mt_tdiv(k,Mnum)=0;
      end
    end
    Mt_flg=0;N_Mt=0;Mt_ivl=0;div_sm=0;
  end
  %% QUESTION: How to set dynamic retention parameter
  % There are two ways to produce daughters
  % (1) Daughter produced by daughters is the 1st generation and have
  % zero damage in it. Daughter produced by mothers becoming
  % daughter will also have zero damage. Mother producing daughters is not
  % the first generation of mother and therefore may contain damage.
    %% daughter division
    Dg_sm_tdiv=sum(sum(ismember(Dg_tdiv,k)));%Daughter total number of tdiv values
    if Dg_sm_tdiv>0
        re1=Altrusm(st);
        re2=re1;
        Dg_flg=1;
        Mt_flg=1;
        for i=1:Dg_sm_tdiv % i is representing Dnum
            %% finding new initial values for 1st time division
            Dg_div_fvl=[DgIn_fvl(k,i) DgDm_fvl(k,i)];
            new_ivl=protein_div(Dg_div_fvl,1);
            Mt_ivl(i,1:2)=new_ivl(1:2);
            Dg_ivl(i,1:2)=new_ivl(3:4);
            nm_Dg_Mt=nm_Dg_Mt+1;% Number of daughters becoming mother
            N_Mt=N_Mt+1;
            N_Dg=N_Dg+1;
        end
        nm_Dgdiv=i;% number of daughters divided
    else
            nm_Dgdiv=0;
    end
    %% Mother division
    Mt_sm_tdiv=sum(sum(ismember(Mt_tdiv,k)));
    if Mt_sm_tdiv>0
        Dg_flg=1;Mt_flg=1;
        for j=1:Mt_sm_tdiv % j is representing Dnum
            if MtDm_fvl(k,j)<Pdeath
                if MtDiv_nm(k,j)>length(re_vec(u,st,:))
                    re1=0;re2=re1;
                else
                    re1=re_vec(u,st,MtDiv_nm(k,j));
                    if re1<0
                        re1=0;
                    end
                    re2=re1;
                end
                MtDiv_nmOld(nm_Dgdiv+j)=MtDiv_nm(k,j);% number of past divisions of particular cell
                Mt_div_fvl=[MtIn_fvl(k,j) MtDm_fvl(k,j)];
                new_ivl=protein_div(Mt_div_fvl,1);
                Mt_ivl(nm_Dgdiv+j,1:2)=new_ivl(1:2);
                Dg_ivl(nm_Dgdiv+j,1:2)=new_ivl(3:4);
                % Now find daughter num of each mother here
                N_Dg=N_Dg+1;
                N_Mt=N_Mt+1;
            elseif nm_Dgdiv>0
                nm_Dgdiv=nm_Dgdiv-1;
            end
        end
    end
    i=0;
end
% Total cells
Dg_alive=ismember(Dg_tdiv,npas:npas+200)+2*ismember(Dg_tdiv,npas-1);
Dg_tot_vert_sm=sum(Dg_alive,2);
TotDg_sm(st,u)=sum(Dg_tot_vert_sm);
Mt_alive=ismember(Mt_tdiv,npas:npas+200)+2*ismember(Mt_tdiv,npas-1);
Mt_tot_vert_sm=sum(Mt_alive,2);
TotMt_sm(st,u)=sum(Mt_tot_vert_sm);
Total_cells(st,u) = TotMt_sm(st,u) + TotDg_sm(st,u);
% for k=1:Total_cells(st,u)
    
% Healthy cells
% Daughter cells damage
valDgDm=DgDm_fvl(dtmax+1:end,:);
sum_valDgDm(st,u)=sum(sum(valDgDm));
% DmDg(st,u)=sum(sum_valDgDm);
% Daughter cells intact
valDgIn=DgIn_fvl(dtmax+1:end,:);
sum_valDgIn(st,u)=sum(sum(valDgIn));
% DmDg(st,u)=sum(sum_valDgIn);

% Mother cells damage
valMtDm=MtDm_fvl(dtmax+1:end,:);
sum_valMtDm(st,u)=sum(sum(valMtDm));
% DmDg(st,u)=sum(sum_valDgDm);
% Mother cells intact
valMtIn=MtIn_fvl(dtmax+1:end,:);
sum_valMtIn(st,u)=sum(sum(valMtIn));

% normMt=MtDm_fvl(dtmax+1:end,:)./Pdeath;
% sort_normMt=sort(normMt(normMt~=0));
% DmMt(st,u)=sum(sort_normMt);
% Ratio function
Tot_Dm(st,u)=sum_valDgDm(st,u)+sum_valMtDm(st,u);
Tot_In(st,u)=sum_valDgIn(st,u)+sum_valMtIn(st,u);
Ratio_DmIn(st,u)=Tot_Dm(st,u)/Tot_In(st,u);
% Ratio from daughter
Tot_DgDm(st,u)=sum_valDgDm(st,u);
Tot_DgIn(st,u)=sum_valDgIn(st,u);
Ratio_DmIn_Dg(st,u)=Tot_DgDm(st,u)/Tot_DgIn(st,u);
% Ratio from mother
Tot_MtDm(st,u)=sum_valMtDm(st,u);
Tot_MtIn(st,u)=sum_valMtIn(st,u);
Ratio_DmIn_Mt(st,u)=Tot_MtDm(st,u)/Tot_MtIn(st,u);
% Payoff daughter
TotDg_payoff(st,u)=Ratio_DmIn_Dg(st,u).*TotDg_sm(st,u);
TotMt_payoff(st,u)=Ratio_DmIn_Mt(st,u).*TotMt_sm(st,u);
TotSizewise_payoff(st,u)=TotDg_payoff(st,u)*Rd+TotMt_payoff(st,u)*Rm;
end
% Make some altruistic effects. If 
if TotSizewise_payoff(1,u)<TotSizewise_payoff(2,u)
    Altrusm(1)=Altrusm(1)*(1-0.1*rand());
    if Altrusm(1)<0
        Altrusm(1)=0;
    end
    if u>1
%     if TotSizewise_payoff(1,u)>TotSizewise_payoff(1,u-1)
%     Altrusm(2)=Altrusm(2)*(1+0.1*rand());
%     if Altrusm(2)>1
%         Altrusm(2)=1;
%     end
    end
else
%     Altrusm(1)=Altrusm(1)*(1+0.1*rand());
%     if Altrusm(1)>1
%         Altrusm(1)=1;
%     end
    Altrusm(2)=Altrusm(2)*(1-0.1*rand());
    if Altrusm(2)<0
        Altrusm(2)=0;
    end
end
vec_Altrusm(:,u)=Altrusm;
end
% check this
% Z=Tot_DgDm./Total_cells;
toc