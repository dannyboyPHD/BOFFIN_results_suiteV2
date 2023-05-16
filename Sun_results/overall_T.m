
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
dir2='/Users/sunbinxuan/Desktop/2E_final_case/2E_t_33100/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   temperature section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename=[exp_direct,'T_3mm.csv'];
exp_T3=load(filename);
filename=[exp_direct,'T_5mm.csv'];
exp_T5=load(filename);
filename=[exp_direct,'T_10mm.csv'];
exp_T10=load(filename);
filename=[exp_direct,'T_20mm.csv'];
exp_T20=load(filename);
filename=[exp_direct,'T_50mm.csv'];
exp_T50=load(filename);
filename=[exp_direct,'T_70mm.csv'];
exp_T70=load(filename);

filename=[dir1,'3mm_T_processed'];
num1_T3=load(filename);
filename=[dir1,'5mm_T_processed'];
num1_T5=load(filename);
filename=[dir1,'10mm_T_processed'];
num1_T10=load(filename);
filename=[dir1,'20mm_T_processed'];
num1_T20=load(filename);
filename=[dir1,'50mm_T_processed'];
num1_T50=load(filename);
filename=[dir1,'70mm_T_processed'];
num1_T70=load(filename);

filename=[dir2,'3mm_T_processed'];
num2_T3=load(filename);
filename=[dir2,'5mm_T_processed'];
num2_T5=load(filename);
filename=[dir2,'10mm_T_processed'];
num2_T10=load(filename);
filename=[dir2,'20mm_T_processed'];
num2_T20=load(filename);
filename=[dir2,'50mm_T_processed'];
num2_T50=load(filename);
filename=[dir2,'70mm_T_processed'];
num2_T70=load(filename);


%filter all the small value to zero
k=size(num1_T3,1);
for ii=1:k
  if(num1_T3(ii,2)<10^-3)
    num1_T3(ii,2)=0; 
  end
  if(num1_T5(ii,2)<10^-3)
    num1_T5(ii,2)=0; 
  end
  if(num1_T10(ii,2)<10^-3)
    num1_T10(ii,2)=0; 
  end  
  if(num1_T20(ii,2)<10^-3)
    num1_T20(ii,2)=0; 
  end
  if(num1_T50(ii,2)<10^-3)
    num1_T50(ii,2)=0; 
  end  
  if(num1_T70(ii,2)<10^-3)
    num1_T70(ii,2)=0; 
  end
  if(num2_T3(ii,2)<10^-3)
    num2_T3(ii,2)=0; 
  end
  if(num2_T5(ii,2)<10^-3)
    num2_T5(ii,2)=0; 
  end
  if(num2_T10(ii,2)<10^-3)
    num2_T10(ii,2)=0; 
  end  
  if(num2_T20(ii,2)<10^-3)
    num2_T20(ii,2)=0; 
  end
  if(num2_T50(ii,2)<10^-3)
    num2_T50(ii,2)=0; 
  end  
  if(num2_T70(ii,2)<10^-3)
    num2_T70(ii,2)=0; 
  end  
 
end

% remove 3/4 num_data to reduce the data density
iii=size(num1_T3,1);
index_remove=zeros(0,1);
for ii=1:iii
  if(mod(ii,2)~=0)
    index_remove=[index_remove,ii];  
  end
end
index_remove(1)=[];
% num1_T3(index_remove,:)=[];
% num1_T5(index_remove,:)=[];
% num1_T10(index_remove,:)=[];
% num1_T20(index_remove,:)=[];
% num1_T50(index_remove,:)=[];
% num1_T70(index_remove,:)=[];
% %
% num2_T3(index_remove,:)=[];
% num2_T5(index_remove,:)=[];
% num2_T10(index_remove,:)=[];
% num2_T20(index_remove,:)=[];
% num2_T50(index_remove,:)=[];
% num2_T70(index_remove,:)=[];
%

figure (1)
width=450;
height=350;
y0=0;
x0=0; 
 p=plot(exp_T3(:,1),exp_T3(:,2),'+k',...
      num1_T3(:,1)*1000,num1_T3(:,2),'-k',... 
      exp_T20(:,1),exp_T20(:,2),'ok',...  
      num1_T20(:,1)*1000,num1_T20(:,2),'--k',... 
      exp_T70(:,1),exp_T70(:,2),'sk',...  
      num1_T70(:,1)*1000,num1_T70(:,2),'-.k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',9);%set(p, {'MarkerFaceColor'} , {'k';'none';'k';'none';'k';'none'})

leg=legend('Exp. HAB 3mm','Num.HAB 3mm','Exp. HAB 20mm','Num.HAB 20mm',...
    'Exp. HAB 70mm','Num.HAB 70mm');set(leg,'FontSize',16,'FontName','Times New Roman','Box','off','Interpreter','latex')

xxx=xlabel('$r$ (mm)');set(xxx,'FontSize',16,'FontName','Times New Roman','Interpreter','latex')
xlim([0,20]);
yticks([300 700 1100 1500 1900 2300]);
ylim([300,2300]);
yyy=ylabel('Temperature (K)');set(yyy,'FontSize',16,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
% savename=[direct_save,'overall_T'];
% saveas(gcf,savename,'epsc');
savename=[direct_save,'overall_T.eps'];
hgexport(gcf,savename);