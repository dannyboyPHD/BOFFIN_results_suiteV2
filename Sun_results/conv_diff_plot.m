% Checking the c/d function
clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Directory
%num_direct='/Users/sunbinxuan/Desktop/debug_cd/10/';
case_num='1';
num_direct='/Users/sunbinxuan/Desktop/workspace_lam/stagbof_lam/les_release/vts/same_TVD_changed_firstpoint/';
x_compen=0.000; %inmersed tube length for compensation

%   Num. et up
Dp_1E=2.0;% primary particle size in 1E method, unit: nm
E=2; %select number of PBE
NPBE=60;% number of pbe elements in input
critical_index=2;
probe_loc=[0.005,0.02];

% Checking the c/d function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Processing num data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Process for property index
num_filename=[num_direct,'central_line.csv'];
num_data = readtable(num_filename,'ReadVariableNames',false);   %read in the data as table
num_data = table2cell(num_data(1,(1:size(num_data,2)-1)));      %trancate theta position
num_data = regexprep(num_data, '\W', '');   %get rid of anything other than character and number
%  extract index number
for i=1:size(num_data,2)
  if(strcmp(num_data(1,i),'Points0'))
    index_x=i;  
  end
  if(strcmp(num_data(1,i),'Points1'))
    index_r=i;  
  end
  if(strcmp(num_data(1,i),'X1'))
    index_fv=i;  
  end
  if(strcmp(num_data(1,i),'drho'))
    index_rho=i;  
  end
  if(strcmp(num_data(1,i),'X0'))
    index_X0=i;  
  end
  if(strcmp(num_data(1,i),'ND0000'))
    index_ND_start=i;  
  end
  if(strcmp(num_data(1,i),'XD0000'))%start of the grid points
    index_XD_start=i;  
  end
  if(strcmp(num_data(1,i),'T'))
    index_T=i;  
  end
end
%  Re-read csv file to generate NUM data storing all the simulation data (except for titles)
num_data = readtable(num_filename,'ReadVariableNames',false);   %read in the data as table
num_data = table2cell(num_data(:,(1:size(num_data,2)-1)));      %transform into cell array and remove theta coordinate
num_data(1,:)=[];                                               %remove the first row (these cells are titles)
m=size(num_data,1);     %number of row
n=size(num_data,2);     %number of column
num_data=str2double(num_data);
%generate PBE grid and prepare other variables
PBE_grid=zeros(NPBE+1,1);
PBE_grid_m=zeros(NPBE,1); %grid mid-point
PBE_dv=zeros(NPBE,1);%interval in volume space
x1_sec=zeros(NPBE,1);%fv for each volume section in volume space
dp_sec=zeros(NPBE,1);%average diameter of primary particles for each volume section in volume space
dp_sec1=zeros(NPBE,1);%test purposes
avNpi_sec=zeros(NPBE,1);%average diameter of primary particles for each volume section in volume space
% for 1E equation
if(E==1)
  for i=1:NPBE
    PBE_grid(i)=num_data(1,index_XD_start+i-1);
    PBE_grid_m(i)=(num_data(1,index_XD_start+i-1)+num_data(1,index_XD_start+i))*0.5;
    PBE_dv(i)=num_data(1,index_XD_start+i)-num_data(1,index_XD_start+i-1);
    if(i==NPBE)
        PBE_grid(i+1)=num_data(1,index_XD_start+i);
    end
  end
  ND_array=zeros(NPBE,1);
  Vp_1E=3.1415926*(Dp_1E*10^-9)^3/6;
  %find critical index
  index_cri_1E=1;
  for i=1:NPBE
      if(Vp_1E<PBE_grid_m(i))
        index_cri_1E=i-1; %after cell i-1. calculate number of primary particles
        break
      end
  end
end
% for 2E equation
if(E==2)
  for i=1:NPBE
    PBE_grid(i)=num_data(1,index_XD_start+i-1);
    PBE_grid_m(i)=(num_data(1,index_XD_start+i-1)+num_data(1,index_XD_start+i))*0.5;
    PBE_dv(i)=num_data(1,index_XD_start+i)-num_data(1,index_XD_start+i-1);
    if(i==NPBE)
        PBE_grid(i+1)=num_data(1,index_XD_start+i);
    end
  end
  ND=index_X0-index_ND_start-1;   %define the total number of pbe scalars
end
diameter_m=(PBE_grid_m*6/3.1415926).^(1/3)*10^9;


%initialise data array
data=zeros(size(probe_loc,2),size(num_data,2));
Ni=zeros(size(probe_loc,2),NPBE);
Npi=zeros(size(probe_loc,2),NPBE);
avn=zeros(size(probe_loc,2),NPBE);
m0=zeros(size(probe_loc,2));
m1=zeros(size(probe_loc,2));
index_probe=zeros(size(probe_loc,2),1);
xloc=num_data(:,index_x);

for kk=1:size(probe_loc,2)
  %locate probe position
  for i=1:m
      if(xloc(i)>probe_loc(kk))
        index_probe(kk)=i-1; %after cell i-1. calculate number of primary particles
        break
      end
  end
  % extract all the info
  data(kk,:)=num_data(index_probe(kk),:);
  %extract number density distribution
  Ni(kk,:)=data(kk,index_ND_start+1:index_ND_start+NPBE);
  %extract primary partcile number density distribution
  Npi(kk,1:critical_index)=data(kk,index_ND_start+1:index_ND_start+critical_index);
  Npi(kk,critical_index+1:NPBE)=data(kk,index_ND_start+NPBE+1:index_ND_start+NPBE+(NPBE-critical_index));
  %calculate avn distribution
  avn(kk,1:critical_index)=1;
  for kkk=critical_index+1:NPBE
    if(Ni(kk,kkk)>0)
      avn(kk,kkk)=Npi(kk,kkk)/Ni(kk,kkk);
    else
      avn(kk,kkk)=0;
    end 
  end
  %calculate moments
  for kkk=1:NPBE
    m0(kk)=m0(kk)+PBE_dv(kkk)*Ni(kk,kkk);
    m1(kk)=m1(kk)+PBE_dv(kkk)*Ni(kk,kkk)*PBE_grid_m(kkk);
  end
end

%plot part
figure (1)
width=450;
height=350;
y0=0;
x0=0; 
%title_plot='Ni distribution';
loglog(PBE_grid_m(:),Ni(1,:),'-k',PBE_grid_m(:),Ni(2,:),'-r',...
       PBE_grid_m(:),Npi(1,:),'-.k',PBE_grid_m(:),Npi(2,:),'-.r',...
     'LineWidth',2,...  
     'MarkerSize',8);
leg=legend('loc1. $n$','loc2. $n$','loc1. $n_p$','loc2. $n_p$');set(leg,'FontSize',20,'FontName','Times New Roman','Box','off','Interpreter','latex')
%t=title(title_plot);set(t,'FontSize',8)
xticks([10^-26 10^-25 10^-24 10^-23 10^-22 10^-21]);
xlim([10^-26,10^-21]);
xxx=xlabel('$v~(\textrm{m}^3) $');set(xxx,'FontSize',20,'FontName','Times New Roman','Interpreter','latex')
yyy=ylabel('Number density $(1/\textrm{m}^{6})$');set(yyy,'FontSize',20,'FontName','Times New Roman','Interpreter','latex')
set(gca,'FontSize',20,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

% %%%%%%%%%%%%%%%
direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';

savename=[direct_save,'test_N_ctr.eps'];
hgexport(gcf,savename);
%saveas(gcf,savename,'epsc');
%hgexport(gcf,savename);
% 
% figure (2)
% width=450;
% height=350;
% y0=0;
% x0=365; 
% title_plot='normalised distribution';
% loglog(PBE_grid_m(:),Ni(1,:)/m0(1),'or',PBE_grid_m(:),Ni(2,:)/m0(2),'-b',...
%      'LineWidth',3, ...
%      'MarkerSize',3);
% leg=legend('loc1.N','loc2.N');set(leg,'FontSize',8)
% t=title(title_plot);set(t,'FontSize',8)
% %xlim([1,90]);
% %ylim([10^-1,100]);
% xl=xlabel('volume ');set(xl,'FontSize',10)
% yl=ylabel('N');set(yl,'FontSize',10)
% set(gca,'FontSize',10)
% set(gcf,'Position',[x0,y0,width,height]);

figure (3)
width=450;
height=350;
y0=0;
x0=365;
%title_plot='avn distribution';
semilogx(PBE_grid_m(:),avn(1,:),'-k',PBE_grid_m(:),avn(2,:),'ok',...
     'LineWidth',2,...  
     'MarkerSize',6);
leg=legend('loc1. $n_{av}$','loc2. $n_{av}$');set(leg,'FontSize',20,'FontName','Times New Roman','Box','off','Interpreter','latex')
%t=title(title_plot);set(t,'FontSize',8)
xticks([10^-26 10^-25 10^-24 10^-23 10^-22 10^-21]);
xlim([10^-26,10^-21]);
yticks([0 2 4 6 8 10 12 14]);
ylim([0,14]);
%ylim([10^-1,100]);
xxx=xlabel('$v~(\textrm{m}^3) $');set(xxx,'FontSize',20,'FontName','Times New Roman','Interpreter','latex')
yyy=ylabel('$n_{av}$');set(yyy,'FontSize',20,'FontName','Times New Roman','Interpreter','latex')
set(gca,'FontSize',20,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
savename=[direct_save,'test_avN_ctr.eps'];
hgexport(gcf,savename);
% saveas(gcf,savename,'epsc');
%savename=[direct_save,'test_avN_ctr.eps'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Process for property index
num_filename=[num_direct,'0_01.csv'];
num_data = readtable(num_filename,'ReadVariableNames',false);   %read in the data as table
num_data = table2cell(num_data(1,(1:size(num_data,2)-1)));      %trancate theta position
num_data = regexprep(num_data, '\W', '');   %get rid of anything other than character and number
%  extract index number
for i=1:size(num_data,2)
  if(strcmp(num_data(1,i),'Points0'))
    index_x=i;  
  end
  if(strcmp(num_data(1,i),'Points1'))
    index_r=i;  
  end
  if(strcmp(num_data(1,i),'X1'))
    index_fv=i;  
  end
  if(strcmp(num_data(1,i),'drho'))
    index_rho=i;  
  end
  if(strcmp(num_data(1,i),'X0'))
    index_X0=i;  
  end
  if(strcmp(num_data(1,i),'ND0000'))
    index_ND_start=i;  
  end
  if(strcmp(num_data(1,i),'XD0000'))%start of the grid points
    index_XD_start=i;  
  end
  if(strcmp(num_data(1,i),'T'))
    index_T=i;  
  end
end
%  Re-read csv file to generate NUM data storing all the simulation data (except for titles)
num_data = readtable(num_filename,'ReadVariableNames',false);   %read in the data as table
num_data = table2cell(num_data(:,(1:size(num_data,2)-1)));      %transform into cell array and remove theta coordinate
num_data(1,:)=[];                                               %remove the first row (these cells are titles)
m=size(num_data,1);     %number of row
n=size(num_data,2);     %number of column
num_data=str2double(num_data);
%generate PBE grid and prepare other variables
PBE_grid=zeros(NPBE+1,1);
PBE_grid_m=zeros(NPBE,1); %grid mid-point
PBE_dv=zeros(NPBE,1);%interval in volume space
x1_sec=zeros(NPBE,1);%fv for each volume section in volume space
dp_sec=zeros(NPBE,1);%average diameter of primary particles for each volume section in volume space
dp_sec1=zeros(NPBE,1);%test purposes
avNpi_sec=zeros(NPBE,1);%average diameter of primary particles for each volume section in volume space
% for 1E equation
if(E==1)
  for i=1:NPBE
    PBE_grid(i)=num_data(1,index_XD_start+i-1);
    PBE_grid_m(i)=(num_data(1,index_XD_start+i-1)+num_data(1,index_XD_start+i))*0.5;
    PBE_dv(i)=num_data(1,index_XD_start+i)-num_data(1,index_XD_start+i-1);
    if(i==NPBE)
        PBE_grid(i+1)=num_data(1,index_XD_start+i);
    end
  end
  ND_array=zeros(NPBE,1);
  Vp_1E=3.1415926*(Dp_1E*10^-9)^3/6;
  %find critical index
  index_cri_1E=1;
  for i=1:NPBE
      if(Vp_1E<PBE_grid_m(i))
        index_cri_1E=i-1; %after cell i-1. calculate number of primary particles
        break
      end
  end
end
% for 2E equation
if(E==2)
  for i=1:NPBE
    PBE_grid(i)=num_data(1,index_XD_start+i-1);
    PBE_grid_m(i)=(num_data(1,index_XD_start+i-1)+num_data(1,index_XD_start+i))*0.5;
    PBE_dv(i)=num_data(1,index_XD_start+i)-num_data(1,index_XD_start+i-1);
    if(i==NPBE)
        PBE_grid(i+1)=num_data(1,index_XD_start+i);
    end
  end
  ND=index_X0-index_ND_start-1;   %define the total number of pbe scalars
end
diameter_m=(PBE_grid_m*6/3.1415926).^(1/3)*10^9;


%initialise data array
data=zeros(size(probe_loc,2),size(num_data,2));
Ni=zeros(size(probe_loc,2),NPBE);
Npi=zeros(size(probe_loc,2),NPBE);
avn=zeros(size(probe_loc,2),NPBE);
m0=zeros(size(probe_loc,2));
m1=zeros(size(probe_loc,2));
index_probe=zeros(size(probe_loc,2),1);
xloc=num_data(:,index_x);

for kk=1:size(probe_loc,2)
  %locate probe position
  for i=1:m
      if(xloc(i)>probe_loc(kk))
        index_probe(kk)=i-1; %after cell i-1. calculate number of primary particles
        break
      end
  end
  % extract all the info
  data(kk,:)=num_data(index_probe(kk),:);
  %extract number density distribution
  Ni(kk,:)=data(kk,index_ND_start+1:index_ND_start+NPBE);
  %extract primary partcile number density distribution
  Npi(kk,1:critical_index)=data(kk,index_ND_start+1:index_ND_start+critical_index);
  Npi(kk,critical_index+1:NPBE)=data(kk,index_ND_start+NPBE+1:index_ND_start+NPBE+(NPBE-critical_index));
  %calculate avn distribution
  avn(kk,1:critical_index)=1;
  for kkk=critical_index+1:NPBE
    if(Ni(kk,kkk)>0)
      avn(kk,kkk)=Npi(kk,kkk)/Ni(kk,kkk);
    else
      avn(kk,kkk)=0;
    end 
  end
  %calculate moments
  for kkk=1:NPBE
    m0(kk)=m0(kk)+PBE_dv(kkk)*Ni(kk,kkk);
    m1(kk)=m1(kk)+PBE_dv(kkk)*Ni(kk,kkk)*PBE_grid_m(kkk);
  end
end
%plot part
% figure (4)
% width=450;
% height=350;
% y0=500;
% x0=0; 
% loglog(PBE_grid_m(:),Ni(1,:),'-k',PBE_grid_m(:),Ni(2,:),'-r',...
%        PBE_grid_m(:),Npi(1,:),'-.k',PBE_grid_m(:),Npi(2,:),'-.r',...
%      'LineWidth',2,...  
%      'MarkerSize',8);
% leg=legend('loc1.N','loc2.N','loc1.Np','loc2.Np');set(leg,'FontSize',20,'FontName','Times New Roman')
% %t=title(title_plot);set(t,'FontSize',8)
% xlim([10^-26,10^-21]);
% xxx=xlabel('V [m^3] ');set(xxx,'FontSize',12,'FontName','Times New Roman')
% yyy=ylabel('N [m^{-6}]');set(yyy,'FontSize',12,'FontName','Times New Roman')
% set(gca,'FontSize',20,'FontName', 'Times New Roman','color','w')
% set(gcf,'Position',[x0,y0,width,height]);
% 
% %plot part
% figure (5)
% width=450;
% height=350;
% y0=500;
% x0=365;
% %title_plot='avn distribution';
% semilogx(PBE_grid_m(:),avn(1,:),'-k',PBE_grid_m(:),avn(2,:),'or',...
%      'LineWidth',2,...  
%      'MarkerSize',6);
% leg=legend('loc1.avN','loc2.avN');set(leg,'FontSize',20,'FontName','Times New Roman')
% %t=title(title_plot);set(t,'FontSize',8)
% xlim([10^-26,10^-21]);
% ylim([0,14]);
% xxx=xlabel('V [m^3] ');set(xxx,'FontSize',12,'FontName','Times New Roman')
% yyy=ylabel('n_{av}');set(yyy,'FontSize',12,'FontName','Times New Roman')
% set(gca,'FontSize',20,'FontName', 'Times New Roman','color','w')
% set(gcf,'Position',[x0,y0,width,height]);
%ti=[case_num,'8.jpg'];
%saveas(gcf,sprintf('%s',ti))
