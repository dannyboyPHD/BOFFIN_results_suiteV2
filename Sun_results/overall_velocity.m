
clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Directory
exp_direct='/Users/sunbinxuan/Desktop/laminar_case/exp_data/santoro_1983/';
axial_location_T=[3,5,10,20,50,70];
axial_location_velocity=[3,5,10,20,40,70,130];
axial_location_fv=[15,50];
axial_location_OH=[7,70];
axial_location_C2H2=[7,20];

%dir1='/Users/sunbinxuan/Desktop/laminar_case/num_data/santoro_1983/BLAN/300K/1E/n6_g5_guo/';
dir1='/Users/sunbinxuan/Desktop/2E_final_case/1E_0.0037/';
%dir2='/Users/sunbinxuan/Desktop/2E_BLANQUART/titania_33100/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   temperature section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=[exp_direct,'U_3mm.csv'];
exp_U3=load(filename);
filename=[exp_direct,'U_5mm.csv'];
exp_U5=load(filename);
filename=[exp_direct,'U_10mm.csv'];
exp_U10=load(filename);
filename=[exp_direct,'U_20mm.csv'];
exp_U20=load(filename);
filename=[exp_direct,'U_40mm.csv'];
exp_U40=load(filename);
filename=[exp_direct,'U_70mm.csv'];
exp_U70=load(filename);

filename=[dir1,'3mm_velocity_processed'];
num1_U3=load(filename);
filename=[dir1,'5mm_velocity_processed'];
num1_U5=load(filename);
filename=[dir1,'10mm_velocity_processed'];
num1_U10=load(filename);
filename=[dir1,'20mm_velocity_processed'];
num1_U20=load(filename);
filename=[dir1,'40mm_velocity_processed'];
num1_U40=load(filename);
filename=[dir1,'70mm_velocity_processed'];
num1_U70=load(filename);

% filename=[dir2,'3mm_velocity_processed'];
% num2_U3=load(filename);
% filename=[dir2,'5mm_velocity_processed'];
% num2_U5=load(filename);
% filename=[dir2,'10mm_velocity_processed'];
% num2_U10=load(filename);
% filename=[dir2,'20mm_velocity_processed'];
% num2_U20=load(filename);
% filename=[dir2,'40mm_velocity_processed'];
% num2_U40=load(filename);
% filename=[dir2,'70mm_velocity_processed'];
% num2_U70=load(filename);


%filter all the small value to zero
k=size(num1_U3,1);
for ii=1:k
  if(num1_U3(ii,2)<10^-3)
    num1_U3(ii,2)=0; 
  end
  if(num1_U5(ii,2)<10^-3)
    num1_U5(ii,2)=0; 
  end
  if(num1_U10(ii,2)<10^-3)
    num1_U10(ii,2)=0; 
  end  
  if(num1_U20(ii,2)<10^-3)
    num1_U20(ii,2)=0; 
  end
  if(num1_U40(ii,2)<10^-3)
    num1_U40(ii,2)=0; 
  end  
  if(num1_U70(ii,2)<10^-3)
    num1_U70(ii,2)=0; 
  end
%   if(num2_U3(ii,2)<10^-3)
%     num2_U3(ii,2)=0; 
%   end
%   if(num2_U5(ii,2)<10^-3)
%     num2_U5(ii,2)=0; 
%   end
%   if(num2_U10(ii,2)<10^-3)
%     num2_U10(ii,2)=0; 
%   end  
%   if(num2_U20(ii,2)<10^-3)
%     num2_U20(ii,2)=0; 
%   end
%   if(num2_U40(ii,2)<10^-3)
%     num2_U40(ii,2)=0; 
%   end  
%   if(num2_U70(ii,2)<10^-3)
%     num2_U70(ii,2)=0; 
%   end  
 
end

% remove 3/4 num_data to reduce the data density
iii=size(num1_U3,1);
index_remove=zeros(0,1);
for ii=1:iii
  if(mod(ii,2)~=0)
    index_remove=[index_remove,ii];  
  end
end
index_remove(1)=[];
num1_U3(index_remove,:)=[];
num1_U5(index_remove,:)=[];
num1_U10(index_remove,:)=[];
num1_U20(index_remove,:)=[];
num1_U40(index_remove,:)=[];
num1_U70(index_remove,:)=[];
%
% num2_U3(index_remove,:)=[];
% num2_U5(index_remove,:)=[];
% num2_U10(index_remove,:)=[];
% num2_U20(index_remove,:)=[];
% num2_U40(index_remove,:)=[];
% num2_U70(index_remove,:)=[];
%

figure (1)
width=450;
height=350;
y0=0;
x0=0; 

 p=plot(exp_U5(:,1),exp_U5(:,2),'dk',...  
      num1_U5(:,1)*1000,num1_U5(:,2),'-k',...
      exp_U20(:,1),exp_U20(:,2),'+k',...  
      num1_U20(:,1)*1000,num1_U20(:,2),'--k',...
      exp_U40(:,1),exp_U40(:,2),'ok',...  
      num1_U40(:,1)*1000,num1_U40(:,2),'-.k',...
      exp_U70(:,1),exp_U70(:,2),'sk',...  
      num1_U70(:,1)*1000,num1_U70(:,2),':k',... 
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',9);%set(p, {'MarkerFaceColor'} , {'k';'none';'k';'none';'k';'none';'k';'none'})

leg=legend('Exp. HAB 5 mm','Num.HAB 5 mm',...
    'Exp. HAB 20 mm','Num.HAB 20 mm',...
    'Exp. HAB 40 mm','Num.HAB 40 mm',...
    'Exp. HAB 70 mm','Num.HAB 70 mm');set(leg,'FontSize',16,'FontName','Times New Roman','Box','off','Interpreter','latex')

%  p=plot(exp_U20(:,1),exp_U20(:,2),'ok',...  
%       num1_U20(:,1)*1000,num1_U20(:,2),'-ok',...
%       num2_U20(:,1)*1000,num2_U20(:,2),'--ok',...     
%       exp_U70(:,1),exp_U70(:,2),'sk',...  
%       num1_U70(:,1)*1000,num1_U70(:,2),'-sk',...
%       num2_U70(:,1)*1000,num2_U70(:,2),'--sk',...  
%      'LineWidth',2,...
%      'MarkerEdgeColor','k',...    
%      'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'k';'none';'none';'k';'none';'none'})
% 
% leg=legend('Exp. data at 20mm','CFV-1PBE at 20mm','CFV-2PBE at 20mm',...
%     'Exp. data at 70mm','CFV-1PBE at 70mm','CFV-2PBE at 70mm');set(leg,'FontSize',16,'FontName','Times New Roman')

xxx=xlabel('$r$ (mm)');set(xxx,'FontSize',16,'FontName','Times New Roman','Interpreter','latex')
xlim([0,20]);
%ylim([0,2500]);
yyy=ylabel('Axial velocity (m/s)');set(yyy,'FontSize',16,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
% savename=[direct_save,'overall_velocity'];
% saveas(gcf,savename,'epsc');
savename=[direct_save,'overall_velocity.eps'];
hgexport(gcf,savename);
