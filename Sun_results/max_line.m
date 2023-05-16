clc
close all
clear all

%deal with maximum soot fv pathline
exp_direct='/Users/sunbinxuan/Desktop/laminar_case/exp_data/santoro_1983/';
num_direct1='/Users/sunbinxuan/Desktop/2E_final_case/1E_0.0037/';
lab1='1E';
num_direct2='/Users/sunbinxuan/Desktop/2E_final_case/2E_s_5nm_100000/';
lab2='2E silica';
num_direct3='/Users/sunbinxuan/Desktop/2E_final_case/2E_t_33100/';
lab3='2E titania';

x_compen=0.000; %inmersed tube length for compensation
weighted_option=3; %1.volume weighted; 2.number density weighted; 3.geometric average
series=1; %selection 1: avn dp, 2:avN, cor dp, 3: cor avn, dp;4: from max_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   process section for max part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process num data
filename=[num_direct1,'max_data'];
num_data1=load(filename);
filename=[num_direct2,'max_data'];
num_data2=load(filename);
filename=[num_direct3,'max_data'];
num_data3=load(filename);

% 1:x loc 2:nuc_fv 3:cond_fv 4:growth_fv 5:oxidation by O2 6:oxidation by OH 
% 7:T 8:area 
filename=[num_direct1,'max_rate'];
num_rate1=load(filename);
filename=[num_direct2,'max_rate'];
num_rate2=load(filename);
filename=[num_direct3,'max_rate'];
num_rate3=load(filename);
% substitute corresponding dp and avNpi based on the input
if(series==1)
  if(weighted_option==1)
    filename=[num_direct1,'max_data_vol'];
    data=load(filename);
    num_data1(:,6)=data(:,1);
    num_data1(:,5)=data(:,2);
    filename=[num_direct2,'max_data_vol'];
    data=load(filename);
    num_data2(:,6)=data(:,1);
    num_data2(:,5)=data(:,2);
    filename=[num_direct3,'max_data_vol'];
    data=load(filename);
    num_data3(:,6)=data(:,1);
    num_data3(:,5)=data(:,2);
  elseif(weighted_option==2)
    filename=[num_direct1,'max_data_num'];
    data=load(filename);
    num_data1(:,6)=data(:,1);
    num_data1(:,5)=data(:,2);
    filename=[num_direct2,'max_data_num'];
    data=load(filename);
    num_data2(:,6)=data(:,1);
    num_data2(:,5)=data(:,2);    
    filename=[num_direct3,'max_data_num'];
    data=load(filename);
    num_data3(:,6)=data(:,1);
    num_data3(:,5)=data(:,2); 
  elseif(weighted_option==3)
    filename=[num_direct1,'max_data_geo'];
    data=load(filename);
    num_data1(:,6)=data(:,1);
    num_data1(:,5)=data(:,2);
    filename=[num_direct2,'max_data_geo'];
    data=load(filename);
    num_data2(:,6)=data(:,1);
    num_data2(:,5)=data(:,2);    
    filename=[num_direct3,'max_data_geo'];
    data=load(filename);
    num_data3(:,6)=data(:,1);
    num_data3(:,5)=data(:,2);  
  end
elseif(series==2)
  if(weighted_option==1)
    filename=[num_direct1,'max_data_vol'];
    data=load(filename);
    num_data1(:,6)=data(:,3);
    num_data1(:,5)=data(:,4);
    filename=[num_direct2,'max_data_vol'];
    data=load(filename);
    num_data2(:,6)=data(:,3);
    num_data2(:,5)=data(:,4);
    filename=[num_direct3,'max_data_vol'];
    data=load(filename);
    num_data3(:,6)=data(:,3);
    num_data3(:,5)=data(:,4);
  elseif(weighted_option==2)
    filename=[num_direct1,'max_data_num'];
    data=load(filename);
    num_data1(:,6)=data(:,3);
    num_data1(:,5)=data(:,4); 
    filename=[num_direct2,'max_data_num'];
    data=load(filename);
    num_data2(:,6)=data(:,3);
    num_data2(:,5)=data(:,4);    
    filename=[num_direct3,'max_data_num'];
    data=load(filename);
    num_data3(:,6)=data(:,3);
    num_data3(:,5)=data(:,4);   
  elseif(weighted_option==3)
    filename=[num_direct1,'max_data_geo'];
    data=load(filename);
    num_data1(:,6)=data(:,3);
    num_data1(:,5)=data(:,4);
    filename=[num_direct2,'max_data_geo'];
    data=load(filename);
    num_data2(:,6)=data(:,3);
    num_data2(:,5)=data(:,4); 
    filename=[num_direct3,'max_data_geo'];
    data=load(filename);
    num_data3(:,6)=data(:,3);
    num_data3(:,5)=data(:,4);  
  end 
elseif(series==3)
  if(weighted_option==1)
    filename=[num_direct1,'max_data_vol'];
    data=load(filename);
    num_data1(:,6)=data(:,5);
    num_data1(:,5)=data(:,6);
    filename=[num_direct2,'max_data_vol'];
    data=load(filename);
    num_data2(:,6)=data(:,5);
    num_data2(:,5)=data(:,6);
    filename=[num_direct3,'max_data_vol'];
    data=load(filename);
    num_data3(:,6)=data(:,5);
    num_data3(:,5)=data(:,6);
  elseif(weighted_option==2)
    filename=[num_direct1,'max_data_num'];
    data=load(filename);
    num_data1(:,6)=data(:,5);
    num_data1(:,5)=data(:,6);
    filename=[num_direct2,'max_data_num'];
    data=load(filename);
    num_data2(:,6)=data(:,5);
    num_data2(:,5)=data(:,6);
    filename=[num_direct3,'max_data_num'];
    data=load(filename);
    num_data3(:,6)=data(:,5);
    num_data3(:,5)=data(:,6);
  elseif(weighted_option==3)
    filename=[num_direct1,'max_data_geo'];
    data=load(filename);
    num_data1(:,6)=data(:,5);
    num_data1(:,5)=data(:,6);
    filename=[num_direct2,'max_data_geo'];
    data=load(filename);
    num_data2(:,6)=data(:,5);
    num_data2(:,5)=data(:,6);    
    filename=[num_direct3,'max_data_geo'];
    data=load(filename);
    num_data3(:,6)=data(:,5);
    num_data3(:,5)=data(:,6);  
  end 
end

% process exp data
filename=[exp_direct,'max_fv'];
exp_fv=load(filename);  %3
filename=[exp_direct,'max_N'];
exp_N=load(filename);   %4
filename=[exp_direct,'max_avNpi'];
exp_np=load(filename);  %5
filename=[exp_direct,'max_dp'];
exp_dp=load(filename);  %6
filename=[exp_direct,'max_Np'];
exp_Nop=load(filename); %7

% remove 1/2 num_data to reduce the data density
iii=size(num_data1,1);
index_remove=zeros(0,1);
for ii=1:iii
  if(mod(ii,6)~=0)
    index_remove=[index_remove,ii];  
  end
end
index_remove(1)=[];
num_data1(index_remove,:)=[];
num_data2(index_remove,:)=[];
num_data3(index_remove,:)=[];

figure (1)
width=450;
height=350;
y0=0;
x0=0; 
p=plot(exp_fv(:,1)*10,exp_fv(:,2),'ok',...
     num_data1(:,1)*1000,num_data1(:,3)*10^6,'-k',...
     num_data2(:,1)*1000,num_data2(:,3)*10^6,'--k',...
     num_data3(:,1)*1000,num_data3(:,3)*10^6,'-.k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'k';'none';'none';'none'})  %none means no filling

leg=legend('Exp.',lab1,lab2,lab3);set(leg,'FontSize',16,'FontName','Times New Roman')

xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([0,90]);
yyy=ylabel('fv [ppm]');set(yyy,'FontSize',12,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

figure (2)
width=450;
height=350;
y0=0;
x0=465; 
p=semilogy(exp_N(:,1),exp_N(:,2)*10^6,'ok',...
     num_data1(:,1)*1000,num_data1(:,4),'-k',...
     num_data2(:,1)*1000,num_data2(:,4),'--k',...
     num_data3(:,1)*1000,num_data3(:,4),'-.k',...     
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'k';'none';'none';'none'})  %none means no filling
leg=legend('Exp.',lab1,lab2,lab3);set(leg,'FontSize',16,'FontName','Times New Roman')
xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([10,80]);
ylim([10^14,10^20]);
yyy=ylabel('N [m^-3]');set(yyy,'FontSize',12,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

figure (3)
width=450;
height=350;
y0=400;
x0=0; 
p=semilogy(exp_np(:,1)*10,exp_np(:,2),'ok',...
     num_data1(:,1)*1000,num_data1(:,6),'-k',...
     num_data2(:,1)*1000,num_data2(:,6),'--k',...
     num_data3(:,1)*1000,num_data3(:,6),'-.k',... 
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'k';'none';'none';'none'})  %none means no filling
leg=legend('Exp.',lab1,lab2,lab3);set(leg,'FontSize',16,'FontName','Times New Roman')
xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([10,60]);
ylim([0,1000]);
yyy=ylabel('np');set(yyy,'FontSize',12,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);


figure (4)
width=450;
height=350;
y0=400;
x0=465; 
p=plot(exp_dp(:,1)*10,exp_dp(:,2),'ok',...
     num_data1(:,1)*1000,num_data1(:,5),'-k',...
     num_data2(:,1)*1000,num_data2(:,5),'--k',...
     num_data3(:,1)*1000,num_data3(:,5),'-.k',...     
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'k';'none';'none';'none'})
leg=legend('Exp.',lab1,lab2,lab3);set(leg,'FontSize',16,'FontName','Times New Roman')
xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([5,80]);
ylim([1,80]);
yyy=ylabel('dp [nm]');set(yyy,'FontSize',12,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

figure (5)
width=450;
height=350;
y0=400;
x0=920; 
p=semilogy(exp_Nop(:,1)*10,exp_Nop(:,2)*10^6,'ok',...
      num_data1(:,1)*1000,num_data1(:,7),'-k',...
      num_data2(:,1)*1000,num_data2(:,7),'--k',...
      num_data3(:,1)*1000,num_data3(:,7),'-.k',...      
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'k';'none';'none';'none'})
leg=legend('Exp.',lab1,lab2,lab3);set(leg,'FontSize',16,'FontName','Times New Roman')
xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([10,60]);
ylim([10^16,10^22]);
yyy=ylabel('Np [m^-3]');set(yyy,'FontSize',12,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

figure (6)
width=450;
height=350;
y0=0;
x0=920;
yyaxis left
p=plot(num_rate1(:,1)*1000,num_rate1(:,4),'-k',...
      num_rate2(:,1)*1000,num_rate2(:,4),'--k',...
      num_rate3(:,1)*1000,num_rate3(:,4),'-.k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'none';'none';'none'})
yl=ylabel('growth rate [m^3/m^3 s]');set(yl,'FontSize',12,'FontName','Times New Roman')
yyaxis right
p=plot(num_rate1(:,1)*1000,num_rate1(:,8),'-r',...
      num_rate2(:,1)*1000,num_rate2(:,8),'--r',...
      num_rate3(:,1)*1000,num_rate3(:,8),'-.r',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'none';'none';'none'})
yr=ylabel('area [m^2/m^3]');set(yr,'FontSize',12,'FontName','Times New Roman')
leg=legend(lab1,lab2,lab3,lab1,lab2,lab3);set(leg,'FontSize',16,'FontName','Times New Roman')
xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([10,60]);
%ylim([10^16,10^22]);
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

figure (7)
width=450;
height=350;
y0=400;
x0=920; 
p=semilogy(num_rate1(:,1)*1000,num_rate1(:,7),'-k',...
      num_rate2(:,1)*1000,num_rate2(:,7),'--k',...
      num_rate3(:,1)*1000,num_rate3(:,7),'-.k',...      
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'none';'none';'none'})
leg=legend(lab1,lab2,lab3);set(leg,'FontSize',16,'FontName','Times New Roman')
xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([10,60]);
%ylim([10^16,10^22]);
yyy=ylabel('T [K]');set(yyy,'FontSize',12,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);