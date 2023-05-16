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
filename=[num_direct1,'integrated_PBE_num'];
num_PBE1=load(filename);
filename=[num_direct2,'integrated_PBE_num'];
num_PBE2=load(filename);
filename=[num_direct3,'integrated_PBE_num'];
num_PBE3=load(filename);
%extract x location
filename=[num_direct1,'integrated_rate_num'];
num_ref=load(filename);
x_loc=num_ref(:,1)*1000;
%extract PBE grid point
filename=[num_direct1,'15mm_PBE_processed'];
num_ref=load(filename);
PBE_dv=num_ref(:,8);   
PBE_point=num_ref(:,1);%61 point, cell boundary
m=size(PBE_dv,1);
PBE_dv(m)=[];          %60 intervels
% get PBE middle point
PBE_m=zeros(m-1,1);    %60 middle point
d_m=zeros(m-1,1);      %60 middle point
for i=1:m-1
  PBE_m(i)=(PBE_point(i)+PBE_point(i+1))*0.5;
  d_m(i)=(PBE_m(i)*6/3.1415926)^(1/3)*10^9;
end
%
for iii=1:size(x_loc,1)
     if((x_loc(iii))>=15)
       index_15mm=iii-1;
       break;
     end
end
for iii=1:size(x_loc,1)
     if((x_loc(iii))>=50)
       index_50mm=iii-1;
       break;
     end
end
for iii=1:size(x_loc,1)
     if((x_loc(iii))>=80)
       index_80mm=iii-1;
       break;
     end
end
%15mm
figure (1)
width=450;
height=350;
y0=0;
x0=0; 
p=loglog(d_m,num_PBE1(index_15mm,:),'-k',...
         d_m,num_PBE2(index_15mm,:),'--k',...
         d_m,num_PBE3(index_15mm,:),'-.k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'},{'none';'none';'none'})  %none means no filling
leg=legend(lab1,lab2,lab3);set(leg,'FontSize',16,'FontName','Times New Roman')
xxx=xlabel('d (nm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([0,120]);
%ylim([0,4]);
yyy=ylabel('integrated number density (1/m^3 m)');set(yyy,'FontSize',12,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

%50mm
figure (2)
width=450;
height=350;
y0=0;
x0=465; 
p=loglog(d_m,num_PBE1(index_50mm,:),'-k',...
         d_m,num_PBE2(index_50mm,:),'--k',...
         d_m,num_PBE3(index_50mm,:),'-.k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'},{'none';'none';'none'})  %none means no filling
leg=legend(lab1,lab2,lab3);set(leg,'FontSize',16,'FontName','Times New Roman')
xxx=xlabel('d (nm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([0,120]);
%ylim([0,4]);
yyy=ylabel('integrated number density (1/m^3 m)');set(yyy,'FontSize',12,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

%80mm
figure (3)
width=450;
height=350;
y0=0;
x0=930; 
p=loglog(d_m,num_PBE1(index_80mm,:),'-k',...
         d_m,num_PBE2(index_80mm,:),'--k',...
         d_m,num_PBE3(index_80mm,:),'-.k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'},{'none';'none';'none'})  %none means no filling
leg=legend(lab1,lab2,lab3);set(leg,'FontSize',16,'FontName','Times New Roman')
xxx=xlabel('d (nm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([0,120]);
%ylim([0,4]);
yyy=ylabel('integrated number density (1/m^3 m)');set(yyy,'FontSize',12,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);
