%% This program finds population on dynamic retention values
clear all
tic
global Rm Rd Pdeath delt npas tmax Pdiv re1   re2 Altrusm remax
tmax=1.5;
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
u=0;
for val=1:12
flag=-1;
v1=0;v2=0;dec=0;
% Altr_st(1:2)=[0.5;0.5];%optim_Altr(iter,:);%
u1=1;u2=1;w1=1;w2=1;
Tot_Dm=0;st1=1;st2=2;
    if val==1
        Altr_st=[0;0];
    elseif val==2
        Altr_st=[1;1];
    else
        figure;u=0;Altr_st=[(val-2)/10;(val-2)/10];
        Tot_Dm_npas=0;DgIn_npas=0;DgDm_npas=0;MtIn_npas=0;MtDm_npas=0;
    end
vec_Altrusm=0;
vec_Altrusm(1:2,1)=Altr_st;
while flag<=5
u=u+1;
for st=st1:st2
re_vec(u,st,1:20)=0;
    if st==1
%         if val==3
            Altrusm=Altr_st(1);
%         end
        Strgy1=ESS_dist(time,rm);%0.875:-0.2:0;%[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];%[0 0.3 0.6 1];%
        div1=Strgy1(end);
        re_vec(u,st,1:div1)=Strgy1(2*div1+1:3*div1);
        re_vec1(w1,1:div1)=re_vec(u,st,1:div1);
        w1=w1+1;
    elseif st==2
%         if val==3
            Altrusm=Altr_st(2);
%         end
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
%     Time = k
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
Tot_Dm_npas(st,u)=(sum(DgDm_npas(st,u,:))+sum(MtDm_npas(st,u,:)));
% Tot_Dm_npas(u)=sum(DgDm_npas(u,:));
% Tot_Dm_npas(u)=sum(MtDm_npas(u,:));
if val>=3
DgDm_len(st,u)=n+length(Dg_ivl(:,1));%72 more
MtDm_len(st,u)=m+length(Mt_ivl(:,1));%72 less
Rat_DgDm_alive(st,u)=sum(DgDm_npas(st,u,:),3)./DgDm_len(st,u);%sum(DgIn_npas(st,u,:))*100;%));DgDm_npas(st,u,:)+
Rat_MtDm_alive(st,u)=sum(MtDm_npas(st,u,:),3)./MtDm_len(st,u);%sum(MtIn_npas(st,u,:))*100;%Pdeath);MtDm_npas(st,u,:)+
Rat_DgMtDm_alive(st,u)=(sum(MtDm_npas(st,u,:))+sum(DgDm_npas(st,u,:)))./Total_cells(st,u);%sum(MtIn_npas(st,u,:))*100;%Pdeath);MtDm_npas(st,u,:)+
Norm_TotDm=(Tot_Dm_npas(st,u)-min_Tot_Dm(st))/(max_Tot_Dm(st)-min_Tot_Dm(st));
Norm_TotCells=1-(Total_cells(st,u)-min_Tot_cells(st))/(max_Tot_cells(st)-min_Tot_cells(st));
if Norm_TotCells<0
    Norm_TotCells
end
if Norm_TotDm<0
    Norm_TotDm
end
Rat_TotDm_alive(st,u)=Norm_TotDm+Norm_TotCells;%Rat_DgMtDm_alive(st,u);%*exp(TotDg_sm(st,u)/Total_cells(st,u))%/100%Rat_MtDm_alive(st,u);%+/+Rat_DgDm_alive(st,u)+
end
end
if val == 1
    max_Tot_Dm(1:2)=221701;%ceil(Tot_Dm_npas(:,u)/10^0)*10^0;%Rat_DgDm_alive(:,u)+Rat_MtDm_alive(:,u);%
    max_Tot_cells(1:2)=1553;%ceil(Total_cells(:,u)/10^0)*10^0;
    break;
elseif val==2
    min_Tot_Dm(1:2)=132961;%floor(Tot_Dm_npas(:,u)/10^0)*10^0;%Rat_DgDm_alive(:,u)+Rat_MtDm_alive(:,u);%
    min_Tot_cells(1:2)=1013;%floor(Total_cells(:,u)/10^0)*10^0;
    break;
end
% (Rat_DgDm_alive(st,u)+Rat_MtDm_alive(st,u))*log(TotDg_sm(st,u))/log(TotMt_sm(st,u))
% log(TotDg_sm(st,u))/log(TotMt_sm(st,u))
% log10(Total_cells(st,u));
if st1>1
    re_vec(u,1,:)=re_vec(u-1,1,:);
    TotDg_sm(1,u)=TotDg_sm(1,u-1);
    TotMt_sm(1,u)=TotMt_sm(1,u-1);
    Total_cells(1,u)=Total_cells(1,u-1);
    MtIn_npas(1,u,:)=MtIn_npas(1,u-1,:);
    MtDm_npas(1,u,:)=MtDm_npas(1,u-1,:);
    DgIn_npas(1,u,:)=DgIn_npas(1,u-1,:);
    DgDm_npas(1,u,:)=DgDm_npas(1,u-1,:);
    DgDm_len(1,u)=DgDm_len(1,u-1);
    MtDm_len(1,u)=MtDm_len(1,u-1);
    Rat_DgDm_alive(1,u)=Rat_DgDm_alive(1,u-1);
    Rat_MtDm_alive(1,u)=Rat_MtDm_alive(1,u-1);
    Rat_TotDm_alive(1,u)=Rat_TotDm_alive(1,u-1);
%     Rat_TotDm_alive(st,u)=Rat_DgDm_alive(st,u)/max(DgDm_len(:,u))*100%Rat_MtDm_alive(st,u);%+/+Rat_DgDm_alive(st,u)+
elseif st2<2
    re_vec(u,2,:)=re_vec(u-1,2,:);
    TotDg_sm(2,u)=TotDg_sm(2,u-1);
    TotMt_sm(2,u)=TotMt_sm(2,u-1);
    Total_cells(2,u)=Total_cells(2,u-1);
    MtIn_npas(2,u,:)=MtIn_npas(2,u-1,:);
    MtDm_npas(2,u,:)=MtDm_npas(2,u-1,:);
    DgIn_npas(2,u,:)=DgIn_npas(2,u-1,:);
    DgDm_npas(2,u,:)=DgDm_npas(2,u-1,:);
    DgDm_len(2,u)=DgDm_len(2,u-1);
    MtDm_len(2,u)=MtDm_len(2,u-1);
    Rat_DgDm_alive(2,u)=Rat_DgDm_alive(2,u-1);
    Rat_MtDm_alive(2,u)=Rat_MtDm_alive(2,u-1);
    Rat_TotDm_alive(2,u)=Rat_TotDm_alive(2,u-1);    
%     Rat_TotDm_alive(st,u)=Rat_DgDm_alive(st,u)/max(DgDm_len(:,u))*100%Rat_MtDm_alive(st,u);%+/+Rat_DgDm_alive(st,u)+
end
% Make some altruistic effects. If
eps=0.02;
Dg_err=eps/DgDm_len(st,u)/20;
Mt_err=eps/MtDm_len(st,u)/20;
err=Dg_err+Mt_err;%+
% err=eps/Total_cells(st,u);
if or(abs(Rat_TotDm_alive(1,u)-Rat_TotDm_alive(2,u))<err,flag>0)
    flag=flag+1
    if and(min(Rat_TotDm_alive(1,:))==Rat_TotDm_alive(1,u-flag+1),min(Rat_TotDm_alive(2,:))==Rat_TotDm_alive(2,u-flag+1))
    if flag<=2
        A_final(flag)=Altr_st(1);
        if flag==1
            Altr_st(1)=A_final(1)+eps*rand();
        elseif flag==2
            if Rat_TotDm_alive(1,u-1)-Rat_TotDm_alive(1,u)>0
                flag=0;
                sprintf('Test Failed! Altrusm(2)= %f has less payoff than %f, set flag = %d',A_final,Altr_st(1),flag)
            else
                Altr_st(1)=A_final(1)-eps*rand();
            end
        end
        if Altr_st(1)>1
            Altr_st(1)=1;
        elseif Altr_st(1)<0
            Altr_st(1)=0;
        end
    elseif flag<=4
        A_final(flag)=Altr_st(2);
        if flag==3
        if Rat_TotDm_alive(1,u-2)-Rat_TotDm_alive(1,u)>0
            flag=0;
            sprintf('Test Failed! Altrusm(2)= %f has less payoff than %f, set flag = %d',A_final,Altr_st(1),flag)
        else
            Altr_st(2)=A_final(3)+eps*rand();
        end
        elseif flag==4
            if Rat_TotDm_alive(2,u-1)-Rat_TotDm_alive(2,u)>0
                flag=0;
                sprintf('Test Failed! Altrusm(2)= %f has less payoff than %f, set flag = %d',A_final,Altr_st(2),flag)
            else
                Altr_st(2)=A_final(3)-eps*rand();
            end
        end
        if Altr_st(2)>1
            Altr_st(2)=1;
        elseif Altr_st(2)<0
            Altr_st(2)=0;
        end
    elseif flag==5
        if Rat_TotDm_alive(2,u-2)-Rat_TotDm_alive(2,u)>0
            Altr_st(1)=A_final(1);
            flag=0;
            sprintf('Test Failed! Altrusm(2)= %f has less payoff than %f, set flag = %d',A_final,Altr_st(1),flag)
        else
            Altr_st(1)=A_final(1);
            Altr_st(2)=A_final(3);
        end
    end
    else
        if not(min(Rat_TotDm_alive(1,:))==Rat_TotDm_alive(1,u))
            Rat_TotDm_alive(1,u)=min(Rat_TotDm_alive(1,:));
            sprintf('Test Failed! Minimum value set strategy 1')
        end
        if not(min(Rat_TotDm_alive(2,:))==Rat_TotDm_alive(2,u))
            Rat_TotDm_alive(2,u)=min(Rat_TotDm_alive(2,:));
            sprintf('Test Failed! Minimum value set strategy 2')
        end
        flag=0;
    end
elseif Rat_TotDm_alive(1,u)>Rat_TotDm_alive(2,u)
    if u==1
        Rat_diff_st1=1;
        Rat_diff_st2=1;
        flag=0;
    else
        Rat_diff_st10=(Rat_TotDm_alive(1,u-1)-Rat_TotDm_alive(1,u));%-ive if prev damage is less
        if Rat_diff_st10<-eps*0.01
            v1=v1+1;%u1=1;
            if mod(v1,2)==1
                Rat_diff_st1=-1;
            else
                Rat_diff_st1=1;
            end
        ratval(1,v1)=Rat_diff_st1;
        end
    end
%     rand_val1=0;
%     while rand_val1<eps/2
        rand_val1=u1*eps*rand();
%     end
    if alt==1
        Altr_st(1)=Altr_st(1)-rand_val1*(Rat_diff_st1/abs(Rat_diff_st1));%
        u1=u1+rand_val1;u2=1;
    else
        Altr_st(1)=Altr_st(1)-rand_val1*(Rat_diff_st1/abs(Rat_diff_st1));%
    end
    if Altr_st(1)>1
        Altr_st(1)=1;
    elseif Altr_st(1)<0
        Altr_st(1)=0;
    end
    st1=1;st2=1;
%     Rat_diff_st1=Rat_TotDm_alive(1,u);
    flag=0;
elseif Rat_TotDm_alive(1,u)<Rat_TotDm_alive(2,u)
    if u==1
        Rat_diff_st1=1;
        Rat_diff_st2=1;
        flag=0;
    else
        Rat_diff_st20=(Rat_TotDm_alive(2,u-1)-Rat_TotDm_alive(2,u));
        if Rat_diff_st20<-eps*0.01
            v2=v2+1;%u2=1;
            if mod(v2,2)==1
                    Rat_diff_st2=-1;
            else
                    Rat_diff_st2=1;
            end
        ratval(2,v2)=Rat_diff_st2;
        end
    end
%     rand_val2=0;
%     while rand_val2<eps/2
        rand_val2=u2*eps*rand();
%     end
    if alt==1
        Altr_st(2)=Altr_st(2)-rand_val2*(Rat_diff_st2/abs(Rat_diff_st2));%
        u2=u2+rand_val2;
        u1=1;
    else
        Altr_st(2)=Altr_st(2)-rand_val2*(Rat_diff_st2/abs(Rat_diff_st2));%
    end
    if Altr_st(2)>1
        Altr_st(2)=1;
    elseif Altr_st(2)<0
        Altr_st(2)=0;
    end
    st1=2;st2=2;
    flag=0;
end
if flag==0
vec_Altrusm(:,dec+u+1)=Altr_st;
plot(vec_Altrusm(1,:),'b-');
hold on;
plot(vec_Altrusm(2,:),'r--');
hold off;
legend('Strategy 1: Distance Measure','Strategy 2: Division Measure','Location','South')
pause(0.001);
else
    dec=dec-1;
end
    if u+dec>5
        if and(sum(abs(diff(vec_Altrusm(1,end-5:end))))==0,sum(abs(diff(vec_Altrusm(2,end-5:end))))==0)
            break;
        end
    end

end
end
% vecAltr_final(iter+1,alt,1:2)=vec_Altrusm(1:2,end);
% vecRatDm_final(iter+1,alt,1:2)=Rat_TotDm_alive(1:2,end);
% end
% [optim_Altr(iter+1,1:2), val]=min(min(vecRatDm_final(iter+1,:,:)));
% iter=iter+1;
% if min(optim_Altr(iter-1,1:2))==min(optim_Altr(iter,1:2))
%     sprintf('Optimum Altrusm values for each strategy is: [S1 S2]=[%f %f]',Altrusm(1),Altrusm(2))
%     break;
% end
% end
% check this
% Z=Tot_DgDm./Total_cells;
toc