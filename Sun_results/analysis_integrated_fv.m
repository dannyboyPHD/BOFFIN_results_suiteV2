clc
close all
clear all

%deal with maximum soot fv pathline
exp_direct='/Users/sunbinxuan/Desktop/laminar_case/exp_data/santoro_1983/';
num_direct1='/Users/sunbinxuan/Desktop/2E_final_case/2E_t_33100/';
lab1='Titania Ea=3.31x10^4';
num_direct2='/Users/sunbinxuan/Desktop/2E_final_case/2E_t_27100/';
lab2='Titania Ea=2.71x10^4';
num_direct3='/Users/sunbinxuan/Desktop/2E_final_case/2E_s_5nm_100000/';
lab3='Silica Ea=10.0x10^4';
num_direct4='/Users/sunbinxuan/Desktop/2E_final_case/2E_s_5nm_80000/';
lab4='Silica Ea=8.0x10^4';

num_direct5='/Users/sunbinxuan/Desktop/2E_final_case/1E_0.0037/';
lab5='1E dc=30.8nm';
num_direct6='/Users/sunbinxuan/Desktop/2E_final_case/1E_0.0037_20nm/';
lab6='1E dc=20.0nm';

x_compen=0.000; %inmersed tube length for compensation
weighted_option=3; %1.volume weighted; 2.number density weighted; 3.geometric average
series=1; %selection 1: avn dp, 2:avN, cor dp, 3: cor avn, dp;4: from max_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   process section for max part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for Integrated fv
filename=[exp_direct,'fv_integrated'];
exp_fv_int=load(filename);
filename=[num_direct1,'integrated_fv_num'];
num_fv1=load(filename);
filename=[num_direct2,'integrated_fv_num'];
num_fv2=load(filename);
filename=[num_direct3,'integrated_fv_num'];
num_fv3=load(filename);
filename=[num_direct4,'integrated_fv_num'];
num_fv4=load(filename);
filename=[num_direct5,'integrated_fv_num'];
num_fv5=load(filename);
filename=[num_direct6,'integrated_fv_num'];
num_fv6=load(filename);
%
figure (1)
width=450;
height=350;
y0=0;
x0=0; 
p=plot(exp_fv_int(:,1)*10,exp_fv_int(:,2)*10^-4*10^10,'ok',...
       num_fv1(:,1)*1000,num_fv1(:,2)*10^10,'-r',...
       num_fv2(:,1)*1000,num_fv2(:,2)*10^10,'--r',...
       num_fv3(:,1)*1000,num_fv3(:,2)*10^10,'-b',...
       num_fv4(:,1)*1000,num_fv4(:,2)*10^10,'--b',...
       num_fv5(:,1)*1000,num_fv5(:,2)*10^10,'-g',...
       num_fv6(:,1)*1000,num_fv6(:,2)*10^10,'--g',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'k';'none';'none';'none';'none';'none';'none'})  %none means no filling

leg=legend('Exp.',lab1,lab2,lab3,lab4,lab5,lab6);set(leg,'FontSize',16,'FontName','Times New Roman')
xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([0,160]);
xticks([0 40 80 120 160]);
ylim([0,2.4]);
yticks([0 0.8 1.6 2.4]);
yyy=ylabel('integrated fv [x10^{-10} m^2]');set(yyy,'FontSize',12,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);