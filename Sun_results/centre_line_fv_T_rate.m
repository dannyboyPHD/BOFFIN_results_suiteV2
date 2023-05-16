clc
close all
clear all

%deal with centre line fv, T, fv_rate and species.

exp_direct='/Users/sunbinxuan/Desktop/laminar_case/exp_data/santoro_1983/';
num_direct1='/Users/sunbinxuan/Desktop/2E_final_case/1E_0.0037/';
lab1='1PBE';
num_direct2='/Users/sunbinxuan/Desktop/2E_final_case/2E_s_5nm_100000/';
lab2='2PBE Silica';
num_direct3='/Users/sunbinxuan/Desktop/2E_final_case/2E_t_33100/';
lab3='2PBE Titania';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Central line T distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%exp
exp_filename=[exp_direct,'central_line_T']; %unit cm, K
exp_data=load(exp_filename);
mmm=size(exp_data,2); %size of the data array
exp_T_x=exp_data(:,1)*10;       %axial location of exp data regarding to T, transform cm to mm
exp_T_y=exp_data(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Central line fv distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_filename=[exp_direct,'central_line_fv']; %unit cm, K
exp_data=load(exp_filename);
mmm=size(exp_data,2); %size of the data array
exp_fv_x=exp_data(:,1)*10;       %axial location of exp data regarding to T, transform cm to mm
exp_fv_y=exp_data(:,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Num.  Central line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_filename=[num_direct1,'central_line_summary']; %unit cm, K
num_data1=load(num_filename); %x,fv,T,avNpi,dp
num_x_1=num_data1(:,1)*1000;
num_fv_1=num_data1(:,2);
num_T_1=num_data1(:,3);
num_avN_1=num_data1(:,4);
num_dp_1=num_data1(:,5);
num_filename=[num_direct2,'central_line_summary']; %unit cm, K
num_data2=load(num_filename); %x,fv,T,avNpi,dp
num_x_2=num_data2(:,1)*1000;
num_fv_2=num_data2(:,2);
num_T_2=num_data2(:,3);
num_avN_2=num_data2(:,4);
num_dp_2=num_data2(:,5);
num_filename=[num_direct3,'central_line_summary']; %unit cm, K
num_data3=load(num_filename); %x,fv,T,avNpi,dp
num_x_3=num_data3(:,1)*1000;
num_fv_3=num_data3(:,2);
num_T_3=num_data3(:,3);
num_avN_3=num_data3(:,4);
num_dp_3=num_data3(:,5);
% rate data
% 1:x location   2:nuc_fv   3: growth_fv 4:nuc  5:pure growth 6:oxidation
% 7:condensation 8:aggregation 9:sintering 10:are number density data at centre line 
% 11: oxidation_fv_O2 12:oxidation_fv_OH 13: condensation_fv
% avNpi and dp along central line
num_filename=[num_direct2,'central_line_rate_summary'];
num_rate=load(num_filename); 
num_rate(:,3)=num_rate(:,5).*num_rate(:,10)./1800;

%species: 1. xloc 2:A2, 3:A2R5, 4:OH, 5:H
num_filename=[num_direct3,'central_line_species'];
num_species=load(num_filename); 


figure (1)
width=450;
height=350;
y0=0;
x0=0; 
p=plot(exp_T_x,exp_T_y,'ok',...
     num_x_1,num_T_1,'-k',...
     num_x_2,num_T_2,'--k',...
     num_x_3,num_T_3,'-.k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'k';'none';'none';'none'})  %none means no filling

leg=legend('Exp.',lab1,lab2,lab3);set(leg,'FontSize',16,'FontName','Times New Roman','Box','off','Interpreter','latex','Location','southeast')

xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',16,'FontName','Times New Roman','Interpreter','latex')
xlim([0,100]);
ylim([200,1800]);
yticks([200 600 1000 1400 1800]);
yyy=ylabel('T (K)');set(yyy,'FontSize',16,'FontName','Times New Roman','Interpreter','latex')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
% savename=[direct_save,'ctr_T'];
% saveas(gcf,savename,'epsc');
savename=[direct_save,'ctr_T.eps'];
hgexport(gcf,savename);

figure (2)
width=450;
height=350;
y0=0;
x0=465; 
p=semilogy(exp_fv_x,exp_fv_y,'ok',...
      num_x_1,num_fv_1,'-k',...
      num_x_2,num_fv_2,'--k',...
      num_x_3,num_fv_3,'-.k',...    
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'k';'none';'none';'none'})  %none means no filling

%leg=legend('Exp.',lab1,lab2,lab3);set(leg,'FontSize',16,'FontName','Times New Roman','Box','off','Interpreter','latex','Location','northwest')

xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',16,'FontName','Times New Roman','Interpreter','latex')
xlim([10,80]);
ylim([0.001,100]);
yyy=ylabel('$f_v$ (ppm)');set(yyy,'FontSize',16,'FontName','Times New Roman','Interpreter','latex')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);
direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
% savename=[direct_save,'ctr_fv'];
% saveas(gcf,savename,'epsc');
savename=[direct_save,'ctr_fv.eps'];
hgexport(gcf,savename);


% remove data to reduce the data density
iii=size(num_rate,1);
index_remove=zeros(0,1);
for ii=1:iii
  if(mod(ii,12)~=0)
    index_remove=[index_remove,ii];  
  end
end
num_rate(index_remove,:)=[];
%

%rate
figure (3)
width=450;
height=350;
y0=400;
x0=0; 
p=plot(num_rate(:,1)*1000,num_rate(:,2),'ok',...
      num_rate(:,1)*1000,num_rate(:,3),'-k',...
      num_rate(:,1)*1000,num_rate(:,13),'xk',...  
      num_rate(:,1)*1000,num_rate(:,12),'-dk',...
      num_rate(:,1)*1000,num_rate(:,11),'-.k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'none';'none';'none';'none';'none'})  %none means no filling
leg=legend('Nucleation','Growth','Condensation','Oxidation by OH','Oxidation by $\textrm{O}_2$');set(leg,'FontSize',15,'FontName','Times New Roman','Box','off','Interpreter','latex')
xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',16,'FontName','Times New Roman','Interpreter','latex')
xlim([10,100]);

yyy=ylabel('$R$ $(\textrm{m}^3/(\textrm{m}^3\cdot s))$');set(yyy,'FontSize',16,'FontName','Times New Roman','Interpreter','latex')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
% savename=[direct_save,'ctr_rate'];
% saveas(gcf,savename,'epsc');
savename=[direct_save,'ctr_rate.eps'];
hgexport(gcf,savename);

%species
% 1. xloc 2:A2, 3:A2R5, 4:OH, 5:H 6:C2H2
% remove 5/6 num_data to reduce the data density
iii=size(num_species,1);
index_remove=zeros(0,1);
for ii=1:iii
  if(mod(ii,10)~=0)
    index_remove=[index_remove,ii];  
  end
end
num_species(index_remove,:)=[];
%
figure (4)
width=450;
height=350;
y0=400;
x0=465; 
p=plot(num_species(:,1)*1000,num_species(:,2)*10,'ok',...
       num_species(:,1)*1000,num_species(:,3)*10,'+k',...
       num_species(:,1)*1000,num_species(:,6)*0.01,'^k',...    
       num_species(:,1)*1000,num_species(:,4),'*k',...  
       num_species(:,1)*1000,num_species(:,5),'dk',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'none';'none';'none';'none';'none'})  %none means no filling
leg=legend('$\textrm{C}_{10}\textrm{H}_8 \times 10$','$\textrm{C}_{12}\textrm{H}_8 \times 10$','$\textrm{C}_2\textrm{H}_2 \times 0.01$','\textrm{OH}','\textrm{H}');set(leg,'FontSize',15,'FontName','Times New Roman','Box','off','Interpreter','latex','Location','northwest')
xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',16,'FontName','Times New Roman','Interpreter','latex')
xlim([0,100]);
ylim([0,0.001]);
yyy=ylabel('Mole fraction');set(yyy,'FontSize',16,'FontName','Times New Roman','Interpreter','latex')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
% savename=[direct_save,'ctr_species'];
% saveas(gcf,savename,'epsc');
savename=[direct_save,'ctr_species.eps'];
hgexport(gcf,savename);
