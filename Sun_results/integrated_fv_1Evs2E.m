% script for soot volume fraction at central line and max fv pathline and
% integrated fv
clc
close all
clear all
%integrated first:

exp_direct='/Users/sunbinxuan/Desktop/laminar_case/exp_data/santoro_1983/';
num_direct1='/Users/sunbinxuan/Desktop/2E_final_case/1E_0.0037/';
lab1='1PBE';
num_direct2='/Users/sunbinxuan/Desktop/2E_final_case/2E_s_5nm_100000/';
lab2='2PBE Silica';
num_direct3='/Users/sunbinxuan/Desktop/2E_final_case/2E_t_33100/';
lab3='2PBE Titania';


% for Integrated fv
filename=[exp_direct,'fv_integrated'];
exp_fv_int=load(filename);
filename=[num_direct1,'integrated_fv_num'];
num_fv1=load(filename);
filename=[num_direct2,'integrated_fv_num'];
num_fv2=load(filename);
filename=[num_direct3,'integrated_fv_num'];
num_fv3=load(filename);
%

figure (1)
width=450;
height=350;
y0=0;
x0=0; 
p=plot(exp_fv_int(:,1)*10,exp_fv_int(:,2)*10^-4*10^10,'ok',...
       num_fv1(:,1)*1000,num_fv1(:,2)*10^10,'-k',...
       num_fv2(:,1)*1000,num_fv2(:,2)*10^10,'--k',...
       num_fv3(:,1)*1000,num_fv3(:,2)*10^10,'-.k',...
     'LineWidth',1,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',9);set(p, {'MarkerFaceColor'} , {'k';'none';'none';'none'})  %none means no filling

leg=legend('Exp.',lab1,lab2,lab3);set(leg,'FontSize',16,'FontName','Times New Roman','Box','off','Interpreter','latex')

xxx=xlabel('HAB (mm)');set(xxx,'FontSize',16,'FontName','Times New Roman','Interpreter','latex')
xlim([0,160]);
xticks([0 40 80 120 160]);
ylim([0,1.8]);
yticks([0 0.6 1.2 1.8]);
yyy=ylabel('Integrated $f_v$ ($\times 10^{-10}~\textrm{m}^2$)');set(yyy,'FontSize',16,'FontName','Times New Roman','Interpreter','latex')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
% savename=[direct_save,'integrated_fv_1Evs2E'];
% saveas(gcf,savename,'epsc');
savename=[direct_save,'integrated_fv_1Evs2E.eps'];
hgexport(gcf,savename);