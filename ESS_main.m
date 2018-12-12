%% Evolutionary stable strategy of Yeast damage rentention mechanism
    %% Decreasing retention by anticipating the possible number of division
    %% Decreasing retention by number of daughters becoming mother.
clear all
global delt npas tmax P RE Rm Rd re1 re2 Altrusm remax div
remax=1;
num_a=10;
a=0:0.1:1;
for j=1:num_a+1
Altrusm(1:1)=a(j);%0.5;
tmax=5;
npas=300*tmax+1;
delt=(tmax/(npas-1));
time=0:delt:tmax;
re1=1;re2=re1;
rm=[0.79];
n=length(rm);
ESS_param();
div=0;
for i=1:n
    Rm=rm(i);Rd=1-Rm;
    divsESS=ESS_divs(time,rm(i));%ESS_divs(time,rm(i));%
    div=divsESS(end);
    parent=divsESS(1:div);
    daugh=divsESS(div+1:2*div);
    RE=divsESS(2*div+1:3*div);
    % Figures
    figure;
    set(gca,'FontSize',20); box on;
    hold on
    plot(1:1:div,parent(1:div),'LineWidth',2);
    plot(1:1:div,daugh(1:div),'LineWidth',2);
    xlabel('Division number','FontSize',20);
    ylabel('Retention re','FontSize',20);
    title('Retention proportion to the mother size','FontSize',20);
    legend1=legend('Retention by mother','Retention by daughter');
    set(legend1,...
    'Position',[0.3571875 0.817803320978987 0.396875 0.047945205479452],...
    'Orientation','horizontal',...
    'FontSize',20);
%    re_{mot}= D_{mot}(t_{div}+1)/D_{tot}(t_{div})','re_{daug}=Rm/Rd*(D_{daug}(t_{div}+1)/D_{tot}(t_{div}))');
    figure
    set(gca,'FontSize',20); box on;
    hold on
    plot(time(1:npas),P(1,1:npas),'LineWidth',2);
    plot(time(1:npas),P(2,1:npas),'LineWidth',2);
    plot(time(1:npas),P(1,1:npas)+P(2,1:npas),'LineWidth',2);
      xlabel('time','FontSize',20);
    ylabel('Intact and Damage Components','FontSize',20);
    legend2=legend('Pint','Pdam','Ptot');
    set(legend2,...
    'Position',[0.3571875 0.817803320978987 0.396875 0.047945205479452],...
    'Orientation','horizontal',...
    'FontSize',20);
end
Div_vec(j)=div;
par_vec(1:div,j)=parent(1:div);
daugh_vec(1:div,j)=daugh(1:div);
Re_vec(1:div,j)=RE(1:div);
par_vec(div+1:end,j)=Rm;
daugh_vec(div+1:end,j)=Rm;
RE=0;div=0;
end
% figure;
% for j=1:num_a+1
%     set(gca,'FontSize',20); box on;
%     hold on
%     plot(1:1:Div_vec(1),par_vec(1:Div_vec(1),j),'LineWidth',2);
%     plot(1:1:Div_vec(1),daugh_vec(1:Div_vec(1),j),'LineWidth',2);
%     xlabel('Division number','FontSize',20);
%     ylabel('Retention re','FontSize',20);
%     title('Retention proportion to the mother size','FontSize',20);
%     legend1=legend('Retention by mother','Retention by daughter');
%     set(legend1,...
%     'Position',[0.3571875 0.817803320978987 0.396875 0.047945205479452],...
%     'Orientation','horizontal',...
%     'FontSize',20);
% %    re_{mot}= D_{mot}(t_{div}+1)/D_{tot}(t_{div})','re_{daug}=Rm/Rd*(D_{daug}(t_{div}+1)/D_{tot}(t_{div}))');
% %     figure
% %     set(gca,'FontSize',20); box on;
% %     hold on
% %     plot(time(1:npas),P(1,1:npas),'LineWidth',2);
% %     plot(time(1:npas),P(2,1:npas),'LineWidth',2);
% %     plot(time(1:npas),P(1,1:npas)+P(2,1:npas),'LineWidth',2);
% %       xlabel('time','FontSize',20);
% %     ylabel('Intact and Damage Components','FontSize',20);
% %     legend2=legend('Pint','Pdam','Ptot');
% %     set(legend2,...
% %     'Position',[0.3571875 0.817803320978987 0.396875 0.047945205479452],...
% %     'Orientation','horizontal',...
% %     'FontSize',20);
% end
% % figure
% % hold on
% % for j=1:num_a+1
% % plot(1:Div_vec(1),Re_vec(:,j),'DisplayName','Re_vec')
% % end