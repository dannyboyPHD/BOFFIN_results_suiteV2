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

axial_location=[15,50,80];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  probe information at three HABs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PBE data array:
% 1: grid point, 2: aggregate number density 3: primary particle number
% density 4:primary particle diamter (nm) 5: average number of primary
% particle per aggregate 6. fv 7: fv percent 8:dv
filename=[num_direct1,'15mm_PBE_processed'];
data=load(filename);
m=size(data,1);
PBE_v1=zeros(size(axial_location,2),m);
PBE_n1=zeros(size(axial_location,2),m);
PBE_np1=zeros(size(axial_location,2),m);
PBE_dp1=zeros(size(axial_location,2),m);
PBE_avn1=zeros(size(axial_location,2),m);
PBE_fv1=zeros(size(axial_location,2),m);
PBE_fv_percent1=zeros(size(axial_location,2),m);
PBE_dv1=zeros(size(axial_location,2),m);
PBE_d1=zeros(size(axial_location,2),m);
PBE_vm1=zeros(size(axial_location,2),m);
PBE_present_x1=zeros(size(axial_location,2),m);
PBE_present_y1=zeros(size(axial_location,2),m);

PBE_v2=zeros(size(axial_location,2),m);
PBE_n2=zeros(size(axial_location,2),m);
PBE_np2=zeros(size(axial_location,2),m);
PBE_dp2=zeros(size(axial_location,2),m);
PBE_avn2=zeros(size(axial_location,2),m);
PBE_fv2=zeros(size(axial_location,2),m);
PBE_fv_percent2=zeros(size(axial_location,2),m);
PBE_dv2=zeros(size(axial_location,2),m);
PBE_d2=zeros(size(axial_location,2),m);
PBE_vm2=zeros(size(axial_location,2),m);
PBE_present_x2=zeros(size(axial_location,2),m);
PBE_present_y2=zeros(size(axial_location,2),m);

PBE_v3=zeros(size(axial_location,2),m);
PBE_n3=zeros(size(axial_location,2),m);
PBE_np3=zeros(size(axial_location,2),m);
PBE_dp3=zeros(size(axial_location,2),m);
PBE_avn3=zeros(size(axial_location,2),m);
PBE_fv3=zeros(size(axial_location,2),m);
PBE_fv_percent3=zeros(size(axial_location,2),m);
PBE_dv3=zeros(size(axial_location,2),m);
PBE_d3=zeros(size(axial_location,2),m);
PBE_vm3=zeros(size(axial_location,2),m);
PBE_present_x3=zeros(size(axial_location,2),m);
PBE_present_y3=zeros(size(axial_location,2),m);

for i=1:size(axial_location,2)
   p_title=num2str(axial_location(i));

   % PBE data array:
   % 1: grid point, 2: aggregate number density 3: primary particle number
   % density 4:primary particle diamter (nm) 5: average number of primary
   % particle per aggregate 6. fv 7: fv percent 8:dv
   filename=[num_direct1,p_title,'mm_PBE_processed'];
   data=load(filename); 
   PBE_v1(i,:)=data(:,1);
   PBE_n1(i,:)=data(:,2);
   PBE_np1(i,:)=data(:,3);
   PBE_dp1(i,:)=data(:,4);
   PBE_avn1(i,:)=data(:,5);
   PBE_fv1(i,:)=data(:,6);
   PBE_fv_percent1(i,:)=data(:,7);
   PBE_dv1(i,:)=data(:,8);
   filename=[num_direct2,p_title,'mm_PBE_processed'];
   data=load(filename); 
   PBE_v2(i,:)=data(:,1);
   PBE_n2(i,:)=data(:,2);
   PBE_np2(i,:)=data(:,3);
   PBE_dp2(i,:)=data(:,4);
   PBE_avn2(i,:)=data(:,5);
   PBE_fv2(i,:)=data(:,6);
   PBE_fv_percent2(i,:)=data(:,7);
   PBE_dv2(i,:)=data(:,8);   
   filename=[num_direct3,p_title,'mm_PBE_processed'];
   data=load(filename); 
   PBE_v3(i,:)=data(:,1);
   PBE_n3(i,:)=data(:,2);
   PBE_np3(i,:)=data(:,3);
   PBE_dp3(i,:)=data(:,4);
   PBE_avn3(i,:)=data(:,5);
   PBE_fv3(i,:)=data(:,6);
   PBE_fv_percent3(i,:)=data(:,7);
   PBE_dv3(i,:)=data(:,8);
   % only pbe n starts from 0:NPBE
   % the rest start from 1:NPBE, 61th = 0.0
   % process it to pbe distribution dN/dlogD
   
   
   %manzolo's way
%    for ii=1:m % m is the number of pbe point, m=61 here
%      PBE_d1(i,ii) = (6/3.141592653*PBE_v1(i,ii))^(1/3); %unit transform from m^3 to m
%      PBE_d2(i,ii) = (6/3.141592653*PBE_v2(i,ii))^(1/3); %unit transform from m^3 to m
%      PBE_d3(i,ii) = (6/3.141592653*PBE_v3(i,ii))^(1/3); %unit transform from m^3 to m
%    end
%    for ii=1:m-1 %m-1 = number of pbe element
%      PBE_present_x1(i,ii)=(PBE_d1(i,ii+1)+PBE_d1(i,ii))*0.5*10^9;   %centre point in diameter space, unit nm
%      PBE_present_x2(i,ii)=(PBE_d2(i,ii+1)+PBE_d2(i,ii))*0.5*10^9;   %centre point in diameter space, unit nm
%      PBE_present_x3(i,ii)=(PBE_d3(i,ii+1)+PBE_d3(i,ii))*0.5*10^9;   %centre point in diameter space, unit nm
%      t=PBE_n1(i,ii+1)*PBE_dv1(i,ii);   %partile number at element centre (per m^3)
%      PBE_present_y1(i,ii) = t/(log10(PBE_d1(i,ii+1))-log10(PBE_d1(i,ii))); %dN/dlog10(D_p)
%      t=PBE_n2(i,ii+1)*PBE_dv2(i,ii);   %partile number at element centre (per m^3)
%      PBE_present_y2(i,ii) = t/(log10(PBE_d2(i,ii+1))-log10(PBE_d2(i,ii))); %dN/dlog10(D_p)
%      t=PBE_n3(i,ii+1)*PBE_dv3(i,ii);   %partile number at element centre (per m^3)
%      PBE_present_y3(i,ii) = t/(log10(PBE_d3(i,ii+1))-log10(PBE_d3(i,ii))); %dN/dlog10(D_p)
%    end
   
%    simple way
   for ii=1:m % m is the number of pbe point, m=61 here
      PBE_d1(i,ii) = (6/3.141592653*PBE_v1(i,ii))^(1/3); %unit transform from m^3 to m
      PBE_d2(i,ii) = (6/3.141592653*PBE_v2(i,ii))^(1/3); %unit transform from m^3 to m
      PBE_d3(i,ii) = (6/3.141592653*PBE_v3(i,ii))^(1/3); %unit transform from m^3 to m
   end
   for ii=1:m-1 %m-1 = number of pbe element
      PBE_present_x1(i,ii)=(PBE_d1(i,ii+1)+PBE_d1(i,ii))*0.5*10^9;   %centre point in diameter space, unit nm
      PBE_present_y1(i,ii) = PBE_n1(i,ii+1);
      PBE_present_x2(i,ii)=(PBE_d2(i,ii+1)+PBE_d2(i,ii))*0.5*10^9;   %centre point in diameter space, unit nm
      PBE_present_y2(i,ii) = PBE_n2(i,ii+1);
      PBE_present_x3(i,ii)=(PBE_d3(i,ii+1)+PBE_d3(i,ii))*0.5*10^9;   %centre point in diameter space, unit nm
      PBE_present_y3(i,ii) = PBE_n3(i,ii+1);
   end
end

figure (1)
width=450;
height=350;
y0=0;
x0=0; 
p=loglog(PBE_present_x1(1,:),PBE_present_y1(1,:),'-k',...
        PBE_present_x2(1,:),PBE_present_y2(1,:),'--k',...
        PBE_present_x3(1,:),PBE_present_y3(1,:),'-.k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'none';'none';'none'})  %none means no filling

leg=legend(lab1,lab2,lab3);set(leg,'FontSize',26,'FontName','Times New Roman','Box','off','Interpreter','latex','Location','southwest')

xxx=xlabel('$d$ (nm)');set(xxx,'FontSize',26,'FontName','Times New Roman','Interpreter','latex')
xlim([1,1000]);
xticks([1 10 100 1000]);
ylim([1e10,1e46]);
yticks([1e10 1e22 1e34 1e46]);
yyy=ylabel('$n(v)~(1/\textrm{m}^6)$');set(yyy,'FontSize',26,'FontName','Times New Roman','Interpreter','latex')
set(gca,'FontSize',26,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);

direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
% savename=[direct_save,'psd_15mm_probe'];
% saveas(gcf,savename,'epsc');
savename=[direct_save,'psd_15mm_probe.eps'];
hgexport(gcf,savename);

figure (2)
width=450;
height=350;
y0=0;
x0=465; 
p=loglog(PBE_present_x1(2,:),PBE_present_y1(2,:),'-k',...
       PBE_present_x2(2,:),PBE_present_y2(2,:),'--k',...
       PBE_present_x3(2,:),PBE_present_y3(2,:),'-.k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'none';'none';'none'})  %none means no filling

%leg=legend(lab1,lab2,lab3);set(leg,'FontSize',26,'FontName','Times New Roman')

xxx=xlabel('$d$ (nm)');set(xxx,'FontSize',26,'FontName','Times New Roman','Interpreter','latex')
xlim([1,1000]);
xticks([1 10 100 1000]);
ylim([1e10,1e46]);
yticks([1e10 1e22 1e34 1e46]);
yyy=ylabel('$n(v)~(1/\textrm{m}^6)$');set(yyy,'FontSize',26,'FontName','Times New Roman','Interpreter','latex')
set(gca,'FontSize',26,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);
direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
% savename=[direct_save,'psd_50mm_probe'];
% saveas(gcf,savename,'epsc');
savename=[direct_save,'psd_50mm_probe.eps'];
hgexport(gcf,savename);

figure (3)
width=450;
height=350;
y0=400;
x0=0; 
p=loglog(PBE_present_x1(3,:),PBE_present_y1(3,:),'-k',...
       PBE_present_x2(3,:),PBE_present_y2(3,:),'--k',...
       PBE_present_x3(3,:),PBE_present_y3(3,:),'-.k',...
     'LineWidth',2,...
     'MarkerEdgeColor','k',...    
     'MarkerSize',8);set(p, {'MarkerFaceColor'} , {'none';'none';'none'})  %none means no filling
%leg=legend(lab1,lab2,lab3);set(leg,'FontSize',26,'FontName','Times New Roman')
xxx=xlabel('$d$ (nm)');set(xxx,'FontSize',26,'FontName','Times New Roman','Interpreter','latex')
xlim([1,1000]);
xticks([1 10 100 1000]);
ylim([1e10,1e46]);
yticks([1e10 1e22 1e34 1e46]);
yyy=ylabel('$n(v)~(1/\textrm{m}^6)$');set(yyy,'FontSize',26,'FontName','Times New Roman','Interpreter','latex')
set(gca,'FontSize',26,'FontName', 'Times New Roman','color','w')
set(gcf,'Position',[x0,y0,width,height]);
direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
% savename=[direct_save,'psd_80mm_probe'];
% saveas(gcf,savename,'epsc');
savename=[direct_save,'psd_80mm_probe.eps'];
hgexport(gcf,savename);