% script for soot volume fraction at central line and max fv pathline and
% integrated fv
clc
close all
clear all
%integrated first:

exp_direct='/Users/sunbinxuan/Desktop/laminar_case/exp_data/santoro_1983/';

num_direct1='/Users/sunbinxuan/Desktop/2E_final_case/1E_0.0037/';
lab1='1E';
num_direct2='/Users/sunbinxuan/Desktop/2E_final_case/2E_s_5nm_100000/';
lab2='2E silica';
num_direct3='/Users/sunbinxuan/Desktop/2E_final_case/2E_t_33100/';
lab3='2E titania';

%1: x loc 2: nuc fv 3:growth fv 4:nuc 5:pure growth kg/m^2 s 6:conden kg/m^2 s 7:oxidation kg/m^2 s
%8:nuc in kg/s.m^3 9: pure growth in kg/s.m^3 10: oxidation in kg/s.m^3
%11: condensation rate in kg/s.m^3 12: area
%13: conden fv m^3/m^3 s 14: oxi fv by O2 15 oxi fv by OH
filename=[num_direct1,'integrated_rate_num'];
num_rate1=load(filename);
filename=[num_direct2,'integrated_rate_num'];
num_rate2=load(filename);
filename=[num_direct3,'integrated_rate_num'];
num_rate3=load(filename);
%
%for sepcies
% 1:x loc 2: C2H2 3: A1-C6H6 4:A2-C10H8 5:A2R5 6:A3-C14H10 7: A4-C16H10
% 8: A4R5  9:P2-C12H10
filename=[num_direct1,'integrated_species_num'];
num_species1=load(filename);
filename=[num_direct2,'integrated_species_num'];
num_species2=load(filename);
filename=[num_direct3,'integrated_species_num'];
num_species3=load(filename);
%

figure (1)
width=450;
height=350;
y0=0;
x0=0; 
p=semilogy(num_species1(:,1)*1000,num_species1(:,4),'--k',...
     num_species1(:,1)*1000,num_species1(:,5),'-.k',...
     num_species1(:,1)*1000,num_species1(:,6),'-.k',...
     num_species1(:,1)*1000,num_species1(:,7),'ok',...
     num_species1(:,1)*1000,num_species1(:,8),'*k',...
     num_species1(:,1)*1000,num_species1(:,9),'^k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'},{'none';'none';'none';'none';'none';'none'})  %none means no filling
leg=legend('A2','A2R5','A3','A4','A4R5','P2');set(leg,'FontSize',16,'FontName','Times New Roman','Box','off','Interpreter','latex')
xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([0,120]);
%ylim([0,4]);
yyy=ylabel('integrated species mole fraction [m^2]');set(yyy,'FontSize',12,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

%
%for sepcies
% 1:x loc 2: C2H2 3: A1-C6H6 4:A2-C10H8 5:A2R5 6:A3-C14H10 7: A4-C16H10
% 8: A4R5  9:P2-C12H10 10:OH 11:H 
figure (2)
width=450;
height=350;
y0=0;
x0=465; 
p=semilogy(num_species1(:,1)*1000,num_species1(:,2),'-k',...
     num_species1(:,1)*1000,num_species1(:,10),'-.k',...
     num_species1(:,1)*1000,num_species1(:,11),'--k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'},{'none';'none';'none'})  %none means no filling
leg=legend('C2H2','OH','H');set(leg,'FontSize',16,'FontName','Times New Roman')
xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([0,120]);
%ylim([0,4]);
yyy=ylabel('integrated species mole fraction [m^2]');set(yyy,'FontSize',12,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);


%1: x loc 2: nuc fv 3:growth fv 4:nuc 5:pure growth kg/m^2 s 6:conden kg/m^2 s 7:oxidation kg/m^2 s
%8:nuc in kg/s.m^3 9: pure growth in kg/s.m^3 10: oxidation in kg/s.m^3
%11: condensation rate in kg/s.m^3 12: area
%13: conden fv m^3/m^3 s 14: oxi fv by O2 15 oxi fv by OH

figure (3)
width=450;
height=350;
y0=0;
x0=930; 
p=plot(num_rate1(:,1)*1000,num_rate1(:,2),'-k',...
     num_rate1(:,1)*1000,num_rate1(:,3),'--k',...
     num_rate1(:,1)*1000,num_rate1(:,13),'-.k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'},{'none';'none';'none'})  %none means no filling
leg=legend('Nucleation','Growth','Condensation');set(leg,'FontSize',16,'FontName','Times New Roman','Box','off','Interpreter','latex')
xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',16,'FontName','Times New Roman','Interpreter','latex')
xlim([0,120]);
ylim([-1.2e-8,1.2e-8]);
yyy=ylabel('Integrated rate $(m^3/(m \cdot s))$');set(yyy,'FontSize',16,'FontName','Times New Roman','Interpreter','latex')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
% savename=[direct_save,'integrated_rate'];
% saveas(gcf,savename,'epsc');
savename=[direct_save,'integrated_rate.eps'];
hgexport(gcf,savename);
% figure (4)
% width=450;
% height=350;
% y0=0;
% x0=930; 
% p=plot(num_rate1(:,1)*1000,num_rate1(:,2),'-k',...
%      num_rate1(:,1)*1000,num_rate1(:,3),'--k',...
%      num_rate1(:,1)*1000,num_rate1(:,13),'-.k',...
%     num_rate2(:,1)*1000,num_rate2(:,2),'-b',...
%      num_rate2(:,1)*1000,num_rate2(:,3),'--b',...
%      num_rate2(:,1)*1000,num_rate2(:,13),'-.b',...
%     num_rate3(:,1)*1000,num_rate3(:,2),'-r',...
%      num_rate3(:,1)*1000,num_rate3(:,3),'--r',...
%      num_rate3(:,1)*1000,num_rate3(:,13),'-.r',...
%      'LineWidth',2,...
%      'MarkerEdgeColor','k',...    
%      'MarkerSize',8);set(p, {'MarkerFaceColor'},{'none';'none';'none';'none';'none';'none';'none';'none';'none'})  %none means no filling
% leg=legend('nuc fv','growth fv','conden fv');set(leg,'FontSize',16,'FontName','Times New Roman')
% xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
% xlim([0,120]);
% %ylim([0,4]);
% yyy=ylabel('integrated rate [m^3/ m s]');set(yyy,'FontSize',12,'FontName','Times New Roman')
% set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
% set(gcf,'Position',[x0,y0,width,height]);
