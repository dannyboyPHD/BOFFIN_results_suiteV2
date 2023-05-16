
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

dir1='/Users/sunbinxuan/Desktop/2E_final_case/1E_0.0037/';
dir2='/Users/sunbinxuan/Desktop/2E_final_case/2E_t_33100/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   temperature section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=[exp_direct,'OH_7mmHAB'];
exp_OH7=load(filename);
filename=[exp_direct,'OH_70mmHAB'];
exp_OH70=load(filename);

filename=[dir1,'7mm_OH_processed'];
num1_OH7=load(filename);
filename=[dir1,'70mm_OH_processed'];
num1_OH70=load(filename);


filename=[dir2,'7mm_OH_processed'];
num2_OH7=load(filename);
filename=[dir2,'70mm_OH_processed'];
num2_OH70=load(filename);


%filter all the small value to zero
k=size(num1_OH7,1);
for ii=1:k
  if(num1_OH7(ii,2)<10^-8)
    num1_OH7(ii,2)=0; 
  end
  if(num1_OH70(ii,2)<10^-8)
    num1_OH70(ii,2)=0; 
  end
  if(num2_OH7(ii,2)<10^-8)
    num2_OH7(ii,2)=0; 
  end
  if(num2_OH70(ii,2)<10^-8)
    num2_OH70(ii,2)=0; 
  end
end  
 
% remove 1/2 num_data to reduce the data density
iii=size(num1_OH7,1);
index_remove=zeros(0,1);
for ii=1:iii
  if(mod(ii,2)~=0)
    index_remove=[index_remove,ii];  
  end
end
index_remove(1)=[];
%num1_OH7(index_remove,:)=[];
%num1_OH70(index_remove,:)=[];
%num2_OH7(index_remove,:)=[];
%num2_OH70(index_remove,:)=[];
%
% remove 1/2 exp_data to reduce the data density
iii=size(exp_OH70,1);
index_remove=zeros(0,1);
for ii=1:iii
  if(mod(ii,3)~=0)
    index_remove=[index_remove,ii];  
  end
end
index_remove(1)=[];
exp_OH70(index_remove,:)=[];

%



figure (1)
width=450;
height=350;
y0=0;
x0=0; 

 p=plot(exp_OH7(:,1)*10,exp_OH7(:,2),'sk',...  
      num1_OH7(:,1)*1000,num1_OH7(:,2),'-k',...
      exp_OH70(:,1)*10,exp_OH70(:,2),'ok',...  
      num1_OH70(:,1)*1000,num1_OH70(:,2),'--k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',9);%set(p, {'MarkerFaceColor'} , {'k';'none';'k';'none'})

leg=legend('Exp. HAB 7mm','Num.HAB 7mm',...
    'Exp. HAB 70mm','Num.HAB 70mm');set(leg,'FontSize',16,'FontName','Times New Roman','Box','off','Interpreter','latex')

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
xlim([0,10]);
ylim([0,0.01]);
yyy=ylabel('OH mole fraction');set(yyy,'FontSize',16,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
% savename=[direct_save,'overall_OH'];
% saveas(gcf,savename,'epsc');
savename=[direct_save,'overall_OH.eps'];
hgexport(gcf,savename);