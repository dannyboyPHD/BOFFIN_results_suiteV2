% for sintering timescale analysis
clc
clear all
close all

x=zeros(100,1);
xaxis=zeros(100,1);
xtaxis=zeros(100,1);
t_s1=zeros(100,1);
t_s2=zeros(100,1);
t_s3=zeros(100,1);
t_s4=zeros(100,1);
t_s5=zeros(100,1);
t_s6=zeros(100,1);
t_s7=zeros(100,1);
t_s8=zeros(100,1);
%black
T1 = 2000;
%blue
T2 = 1600;

%for silica
As_s=1.1*10^-14;
E_s1=10*10^4;
E_s2=11*10^4;

Dc1=2;
Dc1=Dc1*10^-9;
Dc2=5;
Dc2=Dc2*10^-9;
TT1=num2str(T1);
TT2=num2str(T2);

lab1='$E_a=1.0 \times 10^5$ K, $d_c= 2$~nm';
lab2='$E_a=1.1 \times 10^5$ K, $d_c= 2$~nm';
lab3='$E_a=1.0 \times 10^5$ K, $d_c= 5$~nm';


for i=1:100
    xaxis(i)=i;
    x(i)=i*10^-9;
    %for figure 1
    %Tau with different Ea at two T(4 lines), dp=2nm
    t_s1(i)=As_s*x(i)*exp(E_s1/T1*(1-Dc1/x(i)));
    t_s2(i)=As_s*x(i)*exp(E_s2/T1*(1-Dc1/x(i)));
    t_s3(i)=As_s*x(i)*exp(E_s1/T2*(1-Dc1/x(i)));
    t_s4(i)=As_s*x(i)*exp(E_s2/T2*(1-Dc1/x(i)));
    %for figure 2
    %Tau with different dp at two T(4 lines),Ea=10*10^4
    t_s6(i)=As_s*x(i)*exp(E_s1/T1*(1-Dc2/x(i)));
    t_s8(i)=As_s*x(i)*exp(E_s1/T2*(1-Dc2/x(i)));

end
% remove 1/2 num_data to reduce the data density
iii=size(t_s1,1);
index_remove=zeros(0,1);
for ii=7:iii
  if(mod(ii,4)~=0)
    index_remove=[index_remove,ii];  
  end
end
%index_remove(1)=[];
xaxis(index_remove,:)=[];
t_s1(index_remove,:)=[];
t_s2(index_remove,:)=[];
t_s3(index_remove,:)=[];
t_s4(index_remove,:)=[];
t_s6(index_remove,:)=[];
t_s8(index_remove,:)=[];


%for titania
As_t=7.44*10^16;
E_t1=3.31*10^4;
E_t2=2.81*10^4;
t_t1=zeros(100,1);
t_t2=zeros(100,1);
t_t3=zeros(100,1);
t_t4=zeros(100,1);

labt1='$E_a=3.31 \times 10^4$ K';
labt2='$E_a=2.81 \times 10^4$ K';


for i=1:100
    xtaxis(i)=i;
    x(i)=i*10^-9;
    %for figure 2
    %Tau with different Ea at two T(4 lines)
    t_t1(i)=As_t*x(i)^4*T1*exp(E_t1/T1);
    t_t2(i)=As_t*x(i)^4*T1*exp(E_t2/T1);
    t_t3(i)=As_t*x(i)^4*T2*exp(E_t1/T2);
    t_t4(i)=As_t*x(i)^4*T2*exp(E_t2/T2);
end
% remove 1/2 num_data to reduce the data density
iii=size(t_t1,1);
index_remove=zeros(0,1);
for ii=3:iii
  if(mod(ii,4)~=0)
    index_remove=[index_remove,ii];  
  end
end
%index_remove(1)=[];
xtaxis(index_remove,:)=[];
t_t1(index_remove,:)=[];
t_t2(index_remove,:)=[];
t_t3(index_remove,:)=[];
t_t4(index_remove,:)=[];


% figure (1)
% width=450;
% height=350;
% y0=400;
% x0=0;
% loglog(xaxis(:,1),t_s1(:,1),'-+k',xaxis(:,1),t_s2(:,1),'--k',...
%        xaxis(:,1),t_s6(:,1),'-dk',xaxis(:,1),t_s3(:,1),'-+r',...
%        xaxis(:,1),t_s4(:,1),'--r',xaxis(:,1),t_s8(:,1),'-dr',...
%      'LineWidth',2, ...
%      'MarkerSize',7);
% text(40,10e10,'1600K','Color','red','FontSize',16)
% text(50,10e-4,'2000K','Color','black','FontSize',16)
% text(2,2.2e-23,'\leftarrow 2nm','Color','black','FontSize',15)
% text(5,5.5e-23,'\leftarrow 5nm','Color','black','FontSize',15)
% text(1.2,1.0e-2,'\tau = 10^{-4} s','Color','black','FontSize',15)
% line([1,100],[1e-4,1e-4],'Color','black','LineWidth',2,'LineStyle',':')
% %for text arrow
% text(2.5,1.0e2,'Increasing Ea','Color','black','FontSize',15)
% annotation('textarrow',[0.47 0.4],[0.4 0.6])
% text(20,1.0e-14,'Increasing dc','Color','black','FontSize',15)
% annotation('textarrow',[0.55 0.7],[0.65 0.4])
% 
% leg=legend(lab1,lab2,lab3);set(leg,'FontSize',14,'FontName','Times New Roman')
% xxx=xlabel('dp (nm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
% xlim([1,100]);
% %xticks([0.001 0.01 0.1 1 10]);
% ylim([10e-25,10e15]);
% %yticks([10e-25 10e-5 10e10 10e15]);
% 
% yyy=ylabel('\tau [s]');set(yyy,'FontSize',12,'FontName','Times New Roman')
% set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
% set(gcf,'Position',[x0,y0,width,height]);
% 
% 
% figure (2)
% width=450;
% height=350;
% y0=400;
% x0=460;
% loglog(xtaxis(:,1),t_t1(:,1),'-+k',xtaxis(:,1),t_t2(:,1),'--k',...
%        xtaxis(:,1),t_t3(:,1),'-+r',xtaxis(:,1),t_t4(:,1),'--r',...
%      'LineWidth',2, ...
%      'MarkerSize',7);
% text(40,10e0,'1600K','Color','red','FontSize',16)
% text(40,10e-6,'2000K','Color','black','FontSize',16)
% text(1.2,5.0e-4,'\tau = 10^{-4} s','Color','black','FontSize',15)
% line([1,100],[1e-4,1e-4],'Color','black','LineWidth',2,'LineStyle',':')
% text(3.0,1.0e-2,'Increasing Ea','Color','black','FontSize',15)
% annotation('textarrow',[0.6 0.4],[0.3 0.6])
% leg=legend(labt1,labt2);set(leg,'FontSize',15,'FontName','Times New Roman')
% xxx=xlabel('dp (nm) ');set(xxx,'FontSize',12,'FontName','Times New Roman')
% xlim([1,100]);
% ylim([10e-10,10e1]);
% yyy=ylabel('\tau [s]');set(yyy,'FontSize',12,'FontName','Times New Roman')
% set(gca,'FontSize',16,'FontName', 'Times New Roman','color','w')
% set(gcf,'Position',[x0,y0,width,height]);


figure (1)
width=450;
height=350;
y0=400;
x0=0;
loglog(xaxis(:,1),t_s1(:,1),'-+k',xaxis(:,1),t_s2(:,1),'--k',...
       xaxis(:,1),t_s6(:,1),'-dk',xaxis(:,1),t_s3(:,1),'-+r',...
       xaxis(:,1),t_s4(:,1),'--r',xaxis(:,1),t_s8(:,1),'-dr',...
     'LineWidth',2, ...
     'MarkerSize',7);
text(40,10e10,'$1600$ K','Color','red','FontSize',16,'Interpreter','latex')
text(30,10e-8,'$2000$ K','Color','black','FontSize',16,'Interpreter','latex')
text(2.2,2.2e-23,'$\leftarrow$ 2 nm','Color','black','FontSize',15,'Interpreter','latex')
text(5.5,5.5e-23,'$\leftarrow$ 5 nm','Color','black','FontSize',15,'Interpreter','latex')
text(1.2,1.0e-2,'$\tau$ = $10^{-4}$ s','Color','black','FontSize',15,'Interpreter','latex')
line([1,100],[1e-4,1e-4],'Color','black','LineWidth',2,'LineStyle',':')


leg=legend(lab1,lab2,lab3);set(leg,'FontSize',16,'Box','off','Interpreter','latex','FontName','Times New Roman','Location','northwest')
xxx=xlabel('$d_p$ (nm) ');set(xxx,'FontSize',16,'Interpreter','latex','FontName','Times New Roman')
xlim([1,100]);
%xticks([0.001 0.01 0.1 1 10]);
ylim([10e-25,10e15]);
%yticks([10e-25 10e-5 10e10 10e15]);

yyy=ylabel('$\tau$ (s)');set(yyy,'FontSize',16,'Interpreter','latex','FontName','Times New Roman')
set(gca,'FontSize',16,'color','w','FontName','Times New Roman')
set(gcf,'Position',[x0,y0,width,height]);


direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
savename=[direct_save,'Tau_silica'];
saveas(gcf,savename,'epsc');



figure (2)
width=450;
height=350;
y0=400;
x0=460;
loglog(xtaxis(:,1),t_t1(:,1),'-+k',xtaxis(:,1),t_t2(:,1),'--k',...
       xtaxis(:,1),t_t3(:,1),'-+r',xtaxis(:,1),t_t4(:,1),'--r',...
     'LineWidth',2, ...
     'MarkerSize',7);
text(30,10e0,'$1600$ K','Color','red','FontSize',16,'Interpreter','latex')
text(30,10e-6,'$2000$ K','Color','black','FontSize',16,'Interpreter','latex')
text(1.2,5.0e-4,'$\tau = 10^{-4}$ s','Color','black','FontSize',15,'Interpreter','latex')
line([1,100],[1e-4,1e-4],'Color','black','LineWidth',2,'LineStyle',':')

leg=legend(labt1,labt2);set(leg,'FontSize',16,'Box','off','Interpreter','latex','Location','northwest')
xxx=xlabel('$d_p$ (nm) ');set(xxx,'FontSize',16,'Interpreter','latex','FontName','Times New Roman')
xlim([1,100]);
ylim([10e-12,10e1]);
yyy=ylabel('$\tau$ (s)');set(yyy,'FontSize',16,'Interpreter','latex','FontName','Times New Roman')
set(gca,'FontSize',16,'color','w','FontName','Times New Roman')
set(gcf,'Position',[x0,y0,width,height]);

direct_save='/Users/sunbinxuan/Desktop/laminar_case/fig_in_paper/';
savename=[direct_save,'Tau_titania'];
saveas(gcf,savename,'epsc');