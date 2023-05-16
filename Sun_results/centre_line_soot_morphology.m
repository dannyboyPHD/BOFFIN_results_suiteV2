clc
close all
clear all

%deal with centre line fv, T, fv_rate and species.

exp_direct='/Users/sunbinxuan/Desktop/laminar_case/exp_data/santoro_1983/';
num_direct1='/Users/sunbinxuan/Desktop/2E_final_case/1E_0.0037/';
lab1='1E';
num_direct2='/Users/sunbinxuan/Desktop/2E_final_case/2E_s_5nm_100000/';
lab2='2E silica';
num_direct3='/Users/sunbinxuan/Desktop/2E_final_case/2E_t_33100/';
lab3='2E titania';


% select average method for dp and avNpi
% not: 3-2, 11
weighted_option=3; %1.volume weighted; 2.number density weighted; 3.geometric average
series=1; %selection 1: avn dp, 2:avN, cor dp, 3: cor avn, dp;4: from max_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Central line avN distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_filename=[exp_direct,'central_line_avNpi']; %unit cm, K
exp_data=load(exp_filename);
mmm=size(exp_data,2); %size of the data array
exp_avN_x=exp_data(:,1);      
exp_avN_y=exp_data(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Central line dp distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_filename=[exp_direct,'central_line_dp']; %unit cm, K
exp_data=load(exp_filename);
exp_dp_x=exp_data(:,1);      
exp_dp_y=exp_data(:,2);


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

% substitute corresponding dp and avNpi based on the input
% avN, dp; avN,cor dp; cor avN, dp from max_series_data
% 6: average avNpi, 5: average dp in num_data series for ploting
% series selection: 1: avn dp, 2:avN, cor dp, 3: cor avn, dp
if(series==1)
  if(weighted_option==1)
    filename=[num_direct1,'central_line_summary_vol'];
    data=load(filename);
    num_avN_1=data(:,1);
    num_dp_1=data(:,2);
    filename=[num_direct2,'central_line_summary_vol'];
    data=load(filename);
    num_avN_2=data(:,1);
    num_dp_2=data(:,2);
    filename=[num_direct3,'central_line_summary_vol'];
    data=load(filename);
    num_avN_3=data(:,1);
    num_dp_3=data(:,2);
  elseif(weighted_option==2)
    filename=[num_direct1,'central_line_summary_num'];
    data=load(filename);
    num_avN_1=data(:,1);
    num_dp_1=data(:,2);
    filename=[num_direct2,'central_line_summary_num'];
    data=load(filename);
    num_avN_2=data(:,1);
    num_dp_2=data(:,2);
    filename=[num_direct3,'central_line_summary_num'];
    data=load(filename);
    num_avN_3=data(:,1);
    num_dp_3=data(:,2);
  elseif(weighted_option==3)
    filename=[num_direct1,'central_line_summary_geo'];
    data=load(filename);
    num_avN_1=data(:,1);
    num_dp_1=data(:,2);
    filename=[num_direct2,'central_line_summary_geo'];
    data=load(filename);
    num_avN_2=data(:,1);
    num_dp_2=data(:,2);
    filename=[num_direct3,'central_line_summary_geo'];
    data=load(filename);
    num_avN_3=data(:,1);
    num_dp_3=data(:,2);
  end
elseif(series==2)
  if(weighted_option==1)
    filename=[num_direct1,'central_line_summary_vol'];
    data=load(filename);
    num_avN_1=data(:,3);
    num_dp_1=data(:,4);
    filename=[num_direct2,'central_line_summary_vol'];
    data=load(filename);
    num_avN_2=data(:,3);
    num_dp_2=data(:,4);
    filename=[num_direct3,'central_line_summary_vol'];
    data=load(filename);
    num_avN_3=data(:,3);
    num_dp_3=data(:,4);
  elseif(weighted_option==2)
    filename=[num_direct1,'central_line_summary_num'];
    data=load(filename);
    num_avN_1=data(:,3);
    num_dp_1=data(:,4);
    filename=[num_direct2,'central_line_summary_num'];
    data=load(filename);
    num_avN_2=data(:,3);
    num_dp_2=data(:,4);
    filename=[num_direct3,'central_line_summary_num'];
    data=load(filename);
    num_avN_3=data(:,3);
    num_dp_3=data(:,4);
  elseif(weighted_option==3)
    filename=[num_direct1,'central_line_summary_geo'];
    data=load(filename);
    num_avN_1=data(:,3);
    num_dp_1=data(:,4);
    filename=[num_direct2,'central_line_summary_geo'];
    data=load(filename);
    num_avN_2=data(:,3);
    num_dp_2=data(:,4);
    filename=[num_direct3,'central_line_summary_geo'];
    data=load(filename);
    num_avN_3=data(:,3);
    num_dp_3=data(:,4);
  end 
elseif(series==3)    
  if(weighted_option==1)
    filename=[num_direct1,'central_line_summary_vol'];
    data=load(filename);
    num_avN_1=data(:,5);
    num_dp_1=data(:,6);
    filename=[num_direct2,'central_line_summary_vol'];
    data=load(filename);
    num_avN_2=data(:,5);
    num_dp_2=data(:,6);
    filename=[num_direct3,'central_line_summary_vol'];
    data=load(filename);
    num_avN_3=data(:,5);
    num_dp_3=data(:,6);
  elseif(weighted_option==2)
    filename=[num_direct1,'central_line_summary_num'];
    data=load(filename);
    num_avN_1=data(:,5);
    num_dp_1=data(:,6);
    filename=[num_direct2,'central_line_summary_num'];
    data=load(filename);
    num_avN_2=data(:,5);
    num_dp_2=data(:,6);
    filename=[num_direct3,'central_line_summary_num'];
    data=load(filename);
    num_avN_3=data(:,5);
    num_dp_3=data(:,6);
  elseif(weighted_option==3)
    filename=[num_direct1,'central_line_summary_geo'];
    data=load(filename);
    num_avN_1=data(:,5);
    num_dp_1=data(:,6);
    filename=[num_direct2,'central_line_summary_geo'];
    data=load(filename);
    num_avN_2=data(:,5);
    num_dp_2=data(:,6);
    filename=[num_direct3,'central_line_summary_geo'];
    data=load(filename);
    num_avN_3=data(:,5);
    num_dp_3=data(:,6);
  end 
end


figure (1)
width=450;
height=350;
y0=0;
x0=0; 

p=semilogy(exp_avN_x,exp_avN_y,'or',...
           num_x_1,num_avN_1,'-k',...
           num_x_2,num_avN_2,'--k',...
           num_x_3,num_avN_3,'-.k',...    
    'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'none';'none';'none';'none'})  %none means no filling

leg=legend('Exp.',lab1,lab2,lab3);set(leg,'FontSize',14,'FontName','Times New Roman')
xlim([10,90]);
ylim([1,10000]);
xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([0,100]);
yyy=ylabel('number of primary particle per aggregate');set(yyy,'FontSize',12,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

figure (2)
width=450;
height=350;
y0=0;
x0=465; 

p=plot(exp_dp_x,exp_dp_y,'or',...
    num_x_1,num_dp_1,'-k',...
    num_x_2,num_dp_2,'--k',...
    num_x_3,num_dp_3,'-.k',...  
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'none';'none';'none';'none'})  %none means no filling
%leg=legend('Exp.',lab1,lab2,lab3);set(leg,'FontSize',12,'FontName','Times New Roman')
xlim([10,120]);
xxx=xlabel('HAB (mm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
xlim([0,100]);
yyy=ylabel('primary particle diameter (nm)');set(yyy,'FontSize',12,'FontName','Times New Roman')
set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);


