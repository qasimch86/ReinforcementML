%% This program finds population on dynamic retention values
clear all
tic
global Rm Rd Pdeath delt npas tmax Pdiv re1 re2 Altrusm remax
tmax=2;
remax=1;

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
TotSizewise_payoff(:,1)=[10 20];
% while abs(diff(TotSizewise_payoff(:,u)))>.1
ESS_param();
iter=1;alt=1;
vecAltr_final(iter,alt,:)=[0.5;0.5];
optim_Altr(iter,1:2)=min(min(vecAltr_final(iter,:,:)));
% while iter<3
% for alt=1:2
% figure;
u=0;
A0=0; A1=1;dA=.01;
% TotDm_sm_npas=0;DgIn_npas=0;DgDm_npas=0;MtIn_npas=0;MtDm_npas=0;
for Altrusm = A0:dA:A1 %while flag<=5
sprintf('Altruism value: %.3f',Altrusm)
u=u+1;st1=1;st2=2;w1=1;w2=1;
for st=st1:st2
re_vec(u,st,1:20)=0;Mt=0;
    if st==1
%         Altrusm=Altr_st(1);
        Strgy1=ESS_dist(time,rm);%0.875:-0.2:0;%[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];%[0 0.3 0.6 1];%
        div1=Strgy1(end);
        re_vec(u,st,1:div1)=Strgy1(2*div1+1:3*div1);
        re_vec1(w1,1:div1)=re_vec(u,st,1:div1);
        w1=w1+1;
    elseif st==2
%         Altrusm=Altr_st(2);
        Strgy2=ESS_divs(time,rm);
        div2=Strgy2(end);
        re_vec(u,st,1:div2)=Strgy2(2*div2+1:3*div2);
        re_vec2(w2,1:div2)=re_vec(u,st,1:div2);
        w2=w2+1;
    end
IC=[Pdiv*Rd 0];
re_len=length(re_vec(u,st,:));
%% Retention dependence upon anticipated number of possible divisions
    N_Dg=1;N_Mt=0;
    Dg_flg=1;Mt_flg=0;
    nm_Dg_Mt=0;div_sm=0;
    MtDiv_nm=zeros(npas+100,npas);
    Dg_tdiv=0;DgIn_fvl=0;DgDm_fvl=0;
    Mt_tdiv=0;MtIn_fvl=0;MtDm_fvl=0;
    MtDiv_nmOld=0;Mt_ivl=0;Dg_ivl=0;
    Dg_ivl(1,1:2)=IC;nm_Dgdiv=0;
    Dg_n=1;Mt_n=1;MtAlive_ivl=0;DgAlive_ivl=0;
    MtAlive_ivl(1,1:3)=0;DgAlive_ivl(1,1:3)=0;
for k=1:npas-1
     Time = k;
    Mt_ivl_N=0;
    Dg_ivl_N=0;
    %% for daughter
  if Dg_flg==1
    for Dnum=1:N_Dg
        Dg_ivl_N(1:2)=round(Dg_ivl(Dnum,1:2));%Intact and damge comp. of daughter cell 
      if isempty(find(DgIn_givl==Dg_ivl_N(1),1))
          [x0,x1]=min(abs(DgIn_givl-Dg_ivl_N(1)));
          Dg_ivl_srch(1)=find((DgIn_givl==DgIn_givl(x1)));%find intact initial value from general daughter values
          sprintf('Daughter Values changed: from Intact = %.0f to Intact = %.0f',Dg_ivl_N(1),Dg_ivl_srch(1))
      else
          Dg_ivl_srch(1)=find((DgIn_givl==Dg_ivl_N(1)));
      end
      if isempty(find(DgDm_givl==Dg_ivl_N(2),1))
          [y0,y1]=min(abs(DgDm_givl-Dg_ivl_N(2)));
          Dg_ivl_srch(2)=find((DgDm_givl==DgDm_givl(y1)));%find intact initial value from general daughter values
          sprintf('Daughter Values changed: from Damage = %.0f to Damage = %.0f',Dg_ivl_N(2),Dg_ivl_srch(2))
      else
          Dg_ivl_srch(2)=find((DgDm_givl==Dg_ivl_N(2)));
      end
        Dg_tdiv(k,Dnum)=round(scl_fct*Dg_div(Dg_ivl_srch(1),Dg_ivl_srch(2)))+k-1;% Daughter division time
        Dg_tdiv_sm=sum(sum(ismember(Dg_tdiv,Dg_tdiv(k,Dnum))));
        DgIn_fvl(Dg_tdiv(k,Dnum),Dg_tdiv_sm)=DgIn_gfvl(Dg_ivl_srch(1),Dg_ivl_srch(2));% Daughter intact final value
        DgDm_fvl(Dg_tdiv(k,Dnum),Dg_tdiv_sm)=DgDm_gfvl(Dg_ivl_srch(1),Dg_ivl_srch(2));% Daughter damage final value
        if Dg_tdiv(k,Dnum)>=npas
            DgAlive_ivl(Dg_n,1:3)=[k Dg_ivl_N(1) Dg_ivl_N(2)];
            Dg_n=Dg_n+1;
        end
    end
%     DgIn_fvl_vec(k,1:max(Dg_tdiv(k,:)),1:Dg_tdiv_sm)=DgIn_fvl(1:max(Dg_tdiv(k,:)),1:Dg_tdiv_sm);
%     DgDm_fvl_vec(k,1:max(Dg_tdiv(k,:)),1:Dg_tdiv_sm)=DgDm_fvl(1:max(Dg_tdiv(k,:)),1:Dg_tdiv_sm);
    Dg_flg=0;N_Dg=0;Dg_ivl=0;
  end
    %% for mother
  if Mt_flg==1
    for Mnum=1:N_Mt
        Mt_ivl_N(1:2)=round(Mt_ivl(Mnum,1:2));%Intact and damge comp. of daughter cell 
      if isempty(find(MtIn_givl==Mt_ivl_N(1),1))
          [x0,x1]=min(abs(MtIn_givl-Mt_ivl_N(1)));
          Mt_ivl_srch(1)=find((MtIn_givl==MtIn_givl(x1)));
          sprintf('Mother Values changed: from Intact = %.0f to Intact = %.0f',Mt_ivl_N(1),MtIn_givl(x1))
      else
          Mt_ivl_srch(1)=find((MtIn_givl==Mt_ivl_N(1)));
      end
      if isempty(find(MtDm_givl==Mt_ivl_N(2),1))
          [y0,y1]=min(abs(MtDm_givl-Mt_ivl_N(2)));
          Mt_ivl_srch(2)=find((MtDm_givl==MtDm_givl(y1)));
          sprintf('Mother Values changed: from Damage = %.0f to Damage = %.0f',Mt_ivl_N(2),MtDm_givl(y1))
      else
          Mt_ivl_srch(2)=find((MtDm_givl==Mt_ivl_N(2)));
      end
%         Mt_ivl_srch(1)=X;%find intact initial value from general Mother values
%         Mt_ivl_srch(2)=Y;%find intact initial value from general Mother values
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
        if Mt_tdiv(k,Mnum)>=npas
            MtAlive_ivl(Mt_n,1:3)=[k Mt_ivl_N(1) Mt_ivl_N(2)];
            Mt_n=Mt_n+1;
        end
      else
        MtIn_fvl(Mt_tdiv(k,Mnum),Mt_tdiv_sm)=0;% probably when Mt_div = 0
        MtDm_fvl(Mt_tdiv(k,Mnum),Mt_tdiv_sm)=0;% probably when Mt_div = 0
        MtDiv_nm(Mt_tdiv(k,Mnum),Mt_tdiv_sm)=0;
        Mt_tdiv(k,Mnum)=0;
        Mt=Mt+1;
      end
    end
    Mt_flg=0;N_Mt=0;Mt_ivl=0;div_sm=0;
  end
    %% daughter division
    Dg_sm_tdiv=sum(sum(ismember(Dg_tdiv,k)));%Daughter total number of tdiv values
    if Dg_sm_tdiv>0
        re1=re_vec(u,st,1);
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
            else%if nm_Dgdiv>0
                nm_Dgdiv=nm_Dgdiv-1;
            end
        end
        i=0;
    end
    if k>1
    if st==st1
Dg1=ismember(Dg_tdiv,k:k+200)+ismember(Dg_tdiv,k-1)+ismember(Mt_tdiv,k-1);
Dg_sm=sum(Dg1,2);
DgTot_st1(u,k)=sum(Dg_sm);
Mt1=ismember(Mt_tdiv,k:k+200)+ismember(Dg_tdiv,k-1)+ismember(Mt_tdiv,k-1);
Mt_sm1=sum(Mt1,2);
MtTot_st1(u,k)=sum(Mt_sm1);
Pop_st1(u,k)=MtTot_st1(u,k) + DgTot_st1(u,k);
    end
    if st==st2
Dg2=ismember(Dg_tdiv,k:k+200)+ismember(Dg_tdiv,k-1)+ismember(Mt_tdiv,k-1);
Dg_sm2=sum(Dg2,2);
DgTot_st2(u,k)=sum(Dg_sm2);
Mt2=ismember(Mt_tdiv,k:k+200)+ismember(Dg_tdiv,k-1)+ismember(Mt_tdiv,k-1);
Mt_sm2=sum(Mt2,2);
MtTot_st2(u,k)=sum(Mt_sm2);
Pop_st2(u,k)=MtTot_st2(u,k) + DgTot_st2(u,k);
    end
    end
end
% Total cells
Dg_alive=ismember(Dg_tdiv,npas:npas+200)+ismember(Dg_tdiv,npas-1)+ismember(Mt_tdiv,npas-1);
Dg_tot_vert_sm=sum(Dg_alive,2);
TotDg_sm(st,u)=sum(Dg_tot_vert_sm);
Mt_alive=ismember(Mt_tdiv,npas:npas+200)+ismember(Dg_tdiv,npas-1)+ismember(Mt_tdiv,npas-1);
Mt_tot_vert_sm=sum(Mt_alive,2);
TotMt_sm(st,u)=sum(Mt_tot_vert_sm);
Total_cells(st,u) = TotMt_sm(st,u) + TotDg_sm(st,u);
I=0;D=0;
if sum(MtAlive_ivl(1,:))==0
    m=0;
else
for m=1:length(MtAlive_ivl(:,1))
        I(1:2,1)=MtAlive_ivl(m,2:3);
            for k=1:npas-MtAlive_ivl(m,1)
                I(:,k+1)=ESS_desol(I(:,k));
            end
        MtIn_npas(st,u,m)=I(1,k);
        MtDm_npas(st,u,m)=I(2,k);   
end
end
if sum(Mt_ivl(1,:))>0
    MtIn_npas(st,u,m+1:m+length(Mt_ivl(:,1)))=Mt_ivl(:,1);
    MtDm_npas(st,u,m+1:m+length(Mt_ivl(:,1)))=Mt_ivl(:,2);
end
if sum(DgAlive_ivl(1,:))==0
    n=0;
else
for n=1:length(DgAlive_ivl(:,1))
    D(1:2,1)=DgAlive_ivl(n,2:3);
    for k=1:npas-DgAlive_ivl(n,1)
        D(:,k+1)=ESS_desol(D(:,k));
    end
    DgIn_npas(st,u,n)=D(1,k);
    DgDm_npas(st,u,n)=D(2,k);
end
end
if sum(Dg_ivl(1,:))>0
    DgIn_npas(st,u,n+1:n+length(Dg_ivl(:,1)))=Dg_ivl(:,1);
    DgDm_npas(st,u,n+1:n+length(Dg_ivl(:,1)))=Dg_ivl(:,2);
end
DgDm_len(st,u)=n+length(Dg_ivl(:,1));
MtDm_len(st,u)=m+length(Mt_ivl(:,1));
Dm_npas(st,u,1:DgDm_len(st,u))=DgDm_npas(st,u,1:DgDm_len(st,u));
Dm_npas(st,u,DgDm_len(st,u)+1:DgDm_len(st,u)+MtDm_len(st,u))=MtDm_npas(st,u,1:MtDm_len(st,u));
% TotDm_sm_npas(st,u)=sum(DgDm_npas(st,u,:))+sum(MtDm_npas(st,u,:));
end
end
%% minimum and maximum values for Total damage and total cells
TotDm_sm_npas=sum(DgDm_npas,3)+sum(MtDm_npas,3);
max_Tot_Dm=max(max(TotDm_sm_npas/10^0))*10^0;%Rat_DgDm_alive(:,u)+Rat_MtDm_alive(:,u);%
min_Tot_Dm=min(min(TotDm_sm_npas/10^0))*10^0;%Rat_DgDm_alive(:,u)+Rat_MtDm_alive(:,u);%
min_Tot_cells=min(min(Total_cells/10^0))*10^0;
max_Tot_cells=max(max(Total_cells/10^0))*10^0;
% Payoff function
Norm_TotDm=(TotDm_sm_npas-min_Tot_Dm(1))./(max_Tot_Dm(1)-min_Tot_Dm(1));
Norm_TotCells=1-(Total_cells-min_Tot_cells(1))/(max_Tot_cells(1)-min_Tot_cells(1));
Rat_TotDm_alive=(Norm_TotDm)+(Norm_TotCells);%Rat_DgMtDm_alive(st,u);%*exp(TotDg_sm(st,u)/Total_cells(st,u))%/100%Rat_MtDm_alive(st,u);%+/+Rat_DgDm_alive(st,u)+
figure
plot(A0:dA:A1,Rat_TotDm_alive(:,:)');
%% minimum and maximum values for Mother damage and Mother cells
MtDm_sm_npas=sum(MtDm_npas,3);
max_Mt_Dm=max(max(MtDm_sm_npas/10^0))*10^0;%Rat_DgDm_alive(:,u)+Rat_MtDm_alive(:,u);%
max_Mt_cells=max(max(MtDm_len/10^0))*10^0;
min_Mt_Dm=min(min(MtDm_sm_npas/10^0))*10^0;%Rat_DgDm_alive(:,u)+Rat_MtDm_alive(:,u);%
min_Mt_cells=min(min(MtDm_len/10^0))*10^0;
% Payoff function
Norm_MtDm=(MtDm_sm_npas-min_Mt_Dm(1))/(max_Mt_Dm(1)-min_Mt_Dm(1));
Norm_MtCells=1-(MtDm_len-min_Mt_cells(1))/(max_Mt_cells(1)-min_Mt_cells(1));
Rat_MtDm_alive=Norm_MtDm+Norm_MtCells;%Rat_DgMtDm_alive(st,u);%*exp(TotDg_sm(st,u)/Total_cells(st,u))%/100%Rat_MtDm_alive(st,u);%+/+Rat_DgDm_alive(st,u)+
figure
plot(A0:dA:A1,Rat_MtDm_alive(:,:)');
%% minimum and maximu values for Daughter damage and Daughter cells
DgDm_sm_npas=sum(DgDm_npas,3);
max_Dg_Dm=max(max(DgDm_sm_npas/10^0))*10^0;%Rat_DgDm_alive(:,u)+Rat_MtDm_alive(:,u);%
max_Dg_cells=max(max(DgDm_len/10^0))*10^0;
min_Dg_Dm=min(min(DgDm_sm_npas/10^0))*10^0;%Rat_DgDm_alive(:,u)+Rat_MtDm_alive(:,u);%
min_Dg_cells=min(min(DgDm_len/10^0))*10^0;
% Payoff function
Norm_DgDm=(DgDm_sm_npas-min_Dg_Dm(1))/(max_Dg_Dm(1)-min_Dg_Dm(1));
Norm_DgCells=1-(DgDm_len-min_Dg_cells(1))/(max_Dg_cells(1)-min_Dg_cells(1));
Rat_DgDm_alive=Norm_DgCells;%Norm_DgDm;%+Rat_DgMtDm_alive(st,u);%*exp(TotDg_sm(st,u)/Total_cells(st,u))%/100%Rat_MtDm_alive(st,u);%+/+Rat_DgDm_alive(st,u)+
figure
plot(A0:dA:A1,Rat_DgDm_alive(:,:)');
% figure;
% for i=1:1001
% plot(0:delt:tmax-delt,Pop_st1(i,1:Time))
% hold on
% plot(0:delt:tmax-delt,Pop_st2(i,1:Time))
% hold off
% pause(0.01)
% end
%% Damage distribution over mother and daughter cells for each strategy
% damage in Total cells
Dm_npasSt1(:,:)=Dm_npas(1,:,:);
Dm_npasSt2(:,:)=Dm_npas(2,:,:);
maxDm=max(max(max(Dm_npas)));
DL=floor(round(maxDm/100)*100/10);
DA=250;
for A=0:DA:length(A0:dA:A1)
    for k=0:DL:maxDm
        X=and(Dm_npasSt1(A+1,1:Total_cells(1,A+1))>=k,Dm_npasSt1(A+1,1:Total_cells(1,A+1))<k+DL);
        Y=and(Dm_npasSt2(A+1,1:Total_cells(2,A+1))>=k,Dm_npasSt2(A+1,1:Total_cells(2,A+1))<k+DL);
        Dm_dstrSt1(A/DA+1,k/DL+1)=sum(X(1,:))/Total_cells(1,A+1);
        Dm_dstrSt2(A/DA+1,k/DL+1)=sum(Y(1,:))/Total_cells(2,A+1);
    end
X=0;
Y=0;
end
% damage in Mother cells
MtDm_npasSt1(:,:)=MtDm_npas(1,:,:);
MtDm_npasSt2(:,:)=MtDm_npas(2,:,:);
maxMtDm=max(max(max(MtDm_npas)));
DL=floor(round(maxMtDm/100)*100/4);
DA=100;
for A=0:DA:length(A0:dA:A1)
    for k=0:DL:maxMtDm
        X=and(MtDm_npasSt1(A+1,:)>=k,MtDm_npasSt1(A+1,:)<k+DL);
        Y=and(MtDm_npasSt2(A+1,:)>=k,MtDm_npasSt2(A+1,:)<k+DL);
        MtDm_dstrSt1(A/DA+1,k/DL+1)=sum(X(1,:));
        MtDm_dstrSt2(A/DA+1,k/DL+1)=sum(Y(1,:));
    end
end
% damage in daughter cells
DgDm_npasSt1(:,:)=DgDm_npas(1,:,:);
DgDm_npasSt2(:,:)=DgDm_npas(2,:,:);
maxDgDm=max(max(max(DgDm_npas)));
DL=floor(round(maxDgDm/100)*100/4);
DA=100;
for A=0:DA:length(A0:dA:A1)
    for k=0:DL:maxDgDm
        X=and(DgDm_npasSt1(A+1,:)>=k,DgDm_npasSt1(A+1,:)<k+DL);
        Y=and(DgDm_npasSt2(A+1,:)>=k,DgDm_npasSt2(A+1,:)<k+DL);
        DgDm_dstrSt1(A/DA+1,k/DL+1)=sum(X(1,:));
        DgDm_dstrSt2(A/DA+1,k/DL+1)=sum(Y(1,:));
    end
end
Z=0:DL:maxDm;
DgZ=0:DL:maxDm;
MtZ=0:DL:maxDm;
figure
bar(Z,Dm_dstrSt1');
bar(Z,Dm_dstrSt2');
bar(MtZ,MtDm_dstrSt1');
bar(MtZ,MtDm_dstrSt2');
bar(DgZ,DgDm_dstrSt1');
bar(DgZ,DgDm_dstrSt2');
toc