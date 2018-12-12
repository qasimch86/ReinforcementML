%% This program works loads the data file determ_1000_AllPopDistr from the 'data' folder
% The program finds the population distribution by knowing damage, total
% cells, A0, A1 and dA
% For daughter pop distribution, variable names in data file TotDg_sm is replaced by Total_cells
% For daughter pop distribution, variable names in data file DgDm_npas is replaced by Dm_npas
% For Mother pop distribution, variable names in data file TotMt_sm is replaced by Total_cells
% For Mother pop distribution, variable names in data file MtDm_npas is replaced by Dm_npas
clear all
load('data/determ_1000_DgPopDistr.mat');
% load('data/determ_1000_MtPopDistr.mat');
% load('data/determ_1000_AllPopDistr.mat');
dL=10;
dLmax=dL;
clear vars Dm_npasSt1 Dm_npasSt2 DL DA X Y
Dm_npasSt1(:,:)=Dm_npas_t1(1,:,:);
Dm_npasSt2(:,:)=Dm_npas_t1(2,:,:);indi
maxDm_t1=600;%max(max(max(Dm_npas_t1)));
round_maxDm_t1=floor(round(maxDm_t1/100)*100);
DL_t1=round_maxDm_t1/dL;
DA=250;
for A=0:DA:length(A0:dA:A1)
for k=0:DL_t1:round_maxDm_t1
X=and(Dm_npasSt1(A+1,1:Total_cells_t1(1,A+1))>=k,Dm_npasSt1(A+1,1:Total_cells_t1(1,A+1))<k+DL_t1);
Y=and(Dm_npasSt2(A+1,1:Total_cells_t1(2,A+1))>=k,Dm_npasSt2(A+1,1:Total_cells_t1(2,A+1))<k+DL_t1);
Dm_dstrSt1_t1(A/DA+1,k/DL_t1+1)=sum(X(1,:))/Total_cells_t1(1,A+1);
Dm_dstrSt2_t1(A/DA+1,k/DL_t1+1)=sum(Y(1,:))/Total_cells_t1(2,A+1);
end
X=0;
Y=0;
end

clear vars Dm_npasSt1 Dm_npasSt2 DA X Y
Dm_npasSt1(:,:)=Dm_npas_t1p5(1,:,:);
Dm_npasSt2(:,:)=Dm_npas_t1p5(2,:,:);
maxDm_t1p5=600;%max(max(max(Dm_npas_t1p5)));
round_maxDm_t1p5=floor(round(maxDm_t1p5/100)*100);
DL_t1p5=round_maxDm_t1p5/dL;
DA=250;
for A=0:DA:length(A0:dA:A1)
for k=0:DL_t1p5:round_maxDm_t1p5
X=and(Dm_npasSt1(A+1,1:Total_cells_t1p5(1,A+1))>=k,Dm_npasSt1(A+1,1:Total_cells_t1p5(1,A+1))<k+DL_t1p5);
Y=and(Dm_npasSt2(A+1,1:Total_cells_t1p5(2,A+1))>=k,Dm_npasSt2(A+1,1:Total_cells_t1p5(2,A+1))<k+DL_t1p5);
Dm_dstrSt1_t1p5(A/DA+1,k/DL_t1p5+1)=sum(X(1,:))/Total_cells_t1p5(1,A+1);
Dm_dstrSt2_t1p5(A/DA+1,k/DL_t1p5+1)=sum(Y(1,:))/Total_cells_t1p5(2,A+1);
end
X=0;
Y=0;
end

clear vars Dm_npasSt1 Dm_npasSt2 DA X Y
Dm_npasSt1(:,:)=Dm_npas_t1p75(1,:,:);
Dm_npasSt2(:,:)=Dm_npas_t1p75(2,:,:);
maxDm_t1p75=600;%max(max(max(Dm_npas_t1p75)));
round_maxDm_t1p75=floor(round(maxDm_t1p75/100)*100);
DL_t1p75=round_maxDm_t1p75/dL;
DA=250;
for A=0:DA:length(A0:dA:A1)
for k=0:DL_t1p75:round_maxDm_t1p75
X=and(Dm_npasSt1(A+1,1:Total_cells_t1p75(1,A+1))>=k,Dm_npasSt1(A+1,1:Total_cells_t1p75(1,A+1))<k+DL_t1p75);
Y=and(Dm_npasSt2(A+1,1:Total_cells_t1p75(2,A+1))>=k,Dm_npasSt2(A+1,1:Total_cells_t1p75(2,A+1))<k+DL_t1p75);
Dm_dstrSt1_t1p75(A/DA+1,k/DL_t1p75+1)=sum(X(1,:))/Total_cells_t1p75(1,A+1);
Dm_dstrSt2_t1p75(A/DA+1,k/DL_t1p75+1)=sum(Y(1,:))/Total_cells_t1p75(2,A+1);
end
X=0;
Y=0;
end

clear vars Dm_npasSt1 Dm_npasSt2 DA X Y
Dm_npasSt1(:,:)=Dm_npas_t2(1,:,:);
Dm_npasSt2(:,:)=Dm_npas_t2(2,:,:);
maxDm_t2=600;%max(max(max(Dm_npas_t2)));
round_maxDm_t2=floor(round(maxDm_t2/100)*100);
DL_t2=round_maxDm_t2/dL;
DA=250;
for A=0:DA:length(A0:dA:A1)
for k=0:DL_t2:round_maxDm_t2
X=and(Dm_npasSt1(A+1,1:Total_cells_t2(1,A+1))>=k,Dm_npasSt1(A+1,1:Total_cells_t2(1,A+1))<k+DL_t2);
Y=and(Dm_npasSt2(A+1,1:Total_cells_t2(2,A+1))>=k,Dm_npasSt2(A+1,1:Total_cells_t2(2,A+1))<k+DL_t2);
Dm_dstrSt1_t2(A/DA+1,k/DL_t2+1)=sum(X(1,:))/Total_cells_t2(1,A+1);
Dm_dstrSt2_t2(A/DA+1,k/DL_t2+1)=sum(Y(1,:))/Total_cells_t2(2,A+1);
end
X=0;
Y=0;
end

clear vars Dm_npasSt1 Dm_npasSt2 DA X Y
Dm_npasSt1(:,:)=Dm_npas_t2p25(1,:,:);
Dm_npasSt2(:,:)=Dm_npas_t2p25(2,:,:);
maxDm_t2p25=600;%max(max(max(Dm_npas_t2p25)));
round_maxDm_t2p25=floor(round(maxDm_t2p25/100)*100);
DL_t2p25=round_maxDm_t2p25/dL;
DA=250;
for A=0:DA:length(A0:dA:A1)
for k=0:DL_t2p25:round_maxDm_t2p25
X=and(Dm_npasSt1(A+1,1:Total_cells_t2p25(1,A+1))>=k,Dm_npasSt1(A+1,1:Total_cells_t2p25(1,A+1))<k+DL_t2p25);
Y=and(Dm_npasSt2(A+1,1:Total_cells_t2p25(2,A+1))>=k,Dm_npasSt2(A+1,1:Total_cells_t2p25(2,A+1))<k+DL_t2p25);
Dm_dstrSt1_t2p25(A/DA+1,k/DL_t2p25+1)=sum(X(1,:))/Total_cells_t2p25(1,A+1);
Dm_dstrSt2_t2p25(A/DA+1,k/DL_t2p25+1)=sum(Y(1,:))/Total_cells_t2p25(2,A+1);
end
X=0;
Y=0;
end
Z_t1=0:DL_t1:maxDm_t1;
Z_t1p5=0:DL_t1p5:maxDm_t1p5;
Z_t1p75=0:DL_t1p75:maxDm_t1p75;
Z_t2=0:DL_t2:maxDm_t2;
Z_t2p25=0:DL_t2p25:maxDm_t2p25;
clear vars Dm_npasSt1 Dm_npasSt2 DA X Y

% Mean and Standard deviation for time t=1.5, 1.75, 2 and 2.25
cct_St1=cat(3,Dm_dstrSt1_t1p5,Dm_dstrSt1_t1p75,Dm_dstrSt1_t2,Dm_dstrSt1_t2p25);
cct_St2=cat(3,Dm_dstrSt2_t1p5,Dm_dstrSt2_t1p75,Dm_dstrSt2_t2,Dm_dstrSt2_t2p25);
% MeanDm_dstrSt1=(Dm_dstrSt1_t1p5+Dm_dstrSt1_t1p75+Dm_dstrSt1_t2+Dm_dstrSt1_t2p25)/4;
MeanDm_dstrSt1=mean(cct_St1,3);
StdDm_dstrSt1=std(cct_St1,[],3);
MeanDm_dstrSt2=mean(cct_St2,3);
StdDm_dstrSt2=std(cct_St2,[],3);

figure
hold on
bar(Z_t1p75(1:dLmax/dL:dLmax)+30,MeanDm_dstrSt1(:,1:dLmax/dL:dLmax)');
errorbar(Z_t1p75(1:dLmax/dL:dLmax)+12,MeanDm_dstrSt1(1,1:dLmax/dL:dLmax)',StdDm_dstrSt1(1,1:dLmax/dL:dLmax),'.')
errorbar(Z_t1p75(1:dLmax/dL:dLmax)+21,MeanDm_dstrSt1(2,1:dLmax/dL:dLmax)',StdDm_dstrSt1(2,1:dLmax/dL:dLmax),'.')
errorbar(Z_t1p75(1:dLmax/dL:dLmax)+30,MeanDm_dstrSt1(3,1:dLmax/dL:dLmax)',StdDm_dstrSt1(3,1:dLmax/dL:dLmax),'.')
errorbar(Z_t1p75(1:dLmax/dL:dLmax)+39,MeanDm_dstrSt1(4,1:dLmax/dL:dLmax)',StdDm_dstrSt1(4,1:dLmax/dL:dLmax),'.')
errorbar(Z_t1p75(1:dLmax/dL:dLmax)+48,MeanDm_dstrSt1(5,1:dLmax/dL:dLmax)',StdDm_dstrSt1(5,1:dLmax/dL:dLmax),'.')

figure
hold on
bar(Z_t1p75(1:dLmax/dL:dLmax)+30,MeanDm_dstrSt2(:,1:dLmax/dL:dLmax)');
errorbar(Z_t1p75(1:dLmax/dL:dLmax)+12,MeanDm_dstrSt2(1,1:dLmax/dL:dLmax)',StdDm_dstrSt2(1,1:dLmax/dL:dLmax),'.')
errorbar(Z_t1p75(1:dLmax/dL:dLmax)+21,MeanDm_dstrSt2(2,1:dLmax/dL:dLmax)',StdDm_dstrSt2(2,1:dLmax/dL:dLmax),'.')
errorbar(Z_t1p75(1:dLmax/dL:dLmax)+30,MeanDm_dstrSt2(3,1:dLmax/dL:dLmax)',StdDm_dstrSt2(3,1:dLmax/dL:dLmax),'.')
errorbar(Z_t1p75(1:dLmax/dL:dLmax)+39,MeanDm_dstrSt2(4,1:dLmax/dL:dLmax)',StdDm_dstrSt2(4,1:dLmax/dL:dLmax),'.')
errorbar(Z_t1p75(1:dLmax/dL:dLmax)+48,MeanDm_dstrSt2(5,1:dLmax/dL:dLmax)',StdDm_dstrSt2(5,1:dLmax/dL:dLmax),'.')

clear vars Dm_dstrSt1_t1 Dm_dstrSt1_t1p5 Dm_dstrSt1_t1p75 Dm_dstrSt1_t2 Dm_dstrSt1_t2p25
clear vars Dm_dstrSt2_t1 Dm_dstrSt2_t1p5 Dm_dstrSt2_t1p75 Dm_dstrSt2_t2 Dm_dstrSt2_t2p25
