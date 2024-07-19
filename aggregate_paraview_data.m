%% read csv
clear all
%%
test_run_name = 'IS_inf_fusing';
n_steps = '2125000';

cell_data = readtable('celldata_IS_inf_fusing.csv');
grid = readtable('grid_NS_silica.csv');

%% grid size
NX = 200;
NY = 100;
NZ = 1;

total_cells = NX*NY*NZ;

cell_info = zeros(total_cells, 26);
offsets = [1,NX+1, NX+2, (NX+1)*(NY+1),(NX+1)*(NY+1)+1, (NX+1)*(NY+1) + (NX+1), (NX+1)*(NY+1) + (NX+2)];
%% add cellcentre cols to cell_data

cellcentre_0 = zeros(total_cells,1);
cellcentre_1 = zeros(total_cells,1);
cellcentre_2 = zeros(total_cells,1);

cell_data = addvars(cell_data,cellcentre_0);
cell_data = addvars(cell_data,cellcentre_1);
cell_data = addvars(cell_data,cellcentre_2);

%%
% row_no = true cellid +1
cell_data = calc_add_cell_centres(total_cells,offsets,cell_data, cell_info);

% for cell_id = 1:total_cells
% 
%     anchor = [cell_data.StructuredCoordinates_0(cell_id),cell_data.StructuredCoordinates_1(cell_id),cell_data.StructuredCoordinates_2(cell_id)];
%     % find corresponding point data
% 
%     p1 = (grid.StructuredCoordinates_0 == anchor(1)) & (grid.StructuredCoordinates_1 == anchor(2)) & (grid.StructuredCoordinates_2 == anchor(3));
%     p1 = grid(p1,:);
% 
%     anchor_ID = p1.PointID;
% 
%     cell_info(cell_id,1:3) = [p1.Points_0,p1.Points_1,p1.Points_2];
%     
% 
% 
%     for o = 1:length(offsets)
%         p = grid.PointID == anchor_ID + offsets(o);
%         p = grid(p,:);
% 
%         cell_info(cell_id,3*o + 1:3*o + 3) = [p.Points_0, p.Points_1,p.Points_2]; 
%     end
%     
%     % calculate cell centres and add new column to table 
%     % TO DO (1)
%     
%     
%     
%     cell_info(cell_id,25) = cell_info(cell_id,8) - cell_info(cell_id,2);
%     cell_info(cell_id,26) = cell_info(cell_id,8)^2 - cell_info(cell_id,2)^2;
%     
%     [x,y,z] = get_cellcentre(cell_info(cell_id,1:3),cell_info(cell_id,4:6),cell_info(cell_id,7:9),cell_info(cell_id,13:15));
%     cell_data.cellcentre_0(cell_id) = x;
%     cell_data.cellcentre_1(cell_id) = y;
%     cell_data.cellcentre_2(cell_id) = z;
%     
% end

%%
% PSD parameters
m = 60;
index_crit = 13;


% psds = zeros()
ni = zeros(total_cells,m);
np = zeros(total_cells,m);
moments = zeros(total_cells, 4);
rho= zeros(total_cells,1);
dps = zeros(total_cells,2);

[v,dv,v_m, psd_headings] = get_pbe_grid_info(cell_data);


for cell_id = 1:total_cells
    ni(cell_id,:) = get_agg_psd(cell_id,index_crit,m,psd_headings,cell_data);
    np(cell_id,:) = get_np(cell_id,index_crit,m,psd_headings,cell_data);
    rho(cell_id) = get_rho(cell_id,cell_data);
  
    moments(cell_id,1) = sum(  rho(cell_id)*ni(cell_id,:)'.*dv(:)        ); % Mo, zeroth moment
    moments(cell_id,2) = sum(  rho(cell_id)*ni(cell_id,:)'.*dv(:).*v_m(:) ); % M1, first moment
    moments(cell_id,3) = sum(  rho(cell_id)*np(cell_id,:)'.*dv(:)         ); % Mp,o
end

%% calc integrated Soot VF
% get x positions: HAB (m)
HABs = unique(cell_data.cellcentre_0);

int_SVF = calc_ISVF(HABs,cell_data, moments,cell_info);

% int_SVF = zeros(length(HABs),1);
% 
% for h = 1:length(HABs)
%    HAB = HABs(h);
%    
%    rowfilter_at_HAB = cell_data.cellcentre_0 == HAB;
%    cells_at_HAB = cell_data(rowfilter_at_HAB,:);
%    
%    cell_ids = cells_at_HAB.CellID;
%    tmp_SVF = 0;
% %  
%    for c = 1:length(cell_ids)
%       tmp_SVF = tmp_SVF + 0.5*cell_info(cell_ids(c)+1,26)*moments(cell_ids(c)+1,2);
%       
% %       cells_at_HAB.X1(c);cells_at_HAB.X1(c);
%       
% %       cells_at_HAB.X1(c);
% %     ;
%    end
% %    cells_at_HAB.X1(c)
%    int_SVF(h) = tmp_SVF;
%    
% end
%%

exp_intSVF = dlmread('./Sun_results/fv_integrated');

% exp_intSVF = dlmread('./Sun_results/Incipient_Int_fv');
% exp_intSVF = dlmread('./Sun_results/smoking_int_fv');

%%
plot(HABs, 2*pi*int_SVF,'r');
%%
hold on
plot(exp_intSVF(:,1)/100, exp_intSVF(:,2)*10^-4,'o')


%% find annulus of max soot

% [max_soot, maxsoot_cell] = max(moments(:,2)); 
% 
% radius_of_max_soot = cell_data.cellcentre_1(maxsoot_cell);
% % radius_of_max_soot = 0.0026;
%--------------------------------------------
% dp calculation method
%--------------------------------------------
dp_select = 3;
%1 = volume fraction based: M1/Mp,o
%2 = Surface area based:    Surface Area/Mp,o



%3 = 6 Vn/An Diemer approach: 6Vn/An based on point cntact assumption
%------------------------------------------
% extract all cells with that radius

HABs = unique(cell_data.cellcentre_0);

max_soot_cells = zeros(length(HABs),7);


[max_soot_cells,nav_agg_only] = get_path_max_soot(index_crit,HABs,m,ni,np,rho,dv, cell_data,moments,dp_select);


% nav_agg_only = zeros(length(HABs),1);
% 
% for h = 1:length(HABs)
%    HAB = HABs(h);
%    
%    rowfilter_at_HAB = cell_data.cellcentre_0 == HAB;
%    
% %        & cell_data.cellcentre_1 == radius_of_max_soot;
%    
%    cells_at_HAB = cell_data(rowfilter_at_HAB,:);
%    
%    [SVF,max_soot_cell] = max(cells_at_HAB.X1);
%    
%    max_soot_cells(h,1) = cells_at_HAB.CellID(max_soot_cell);
%    max_soot_cells(h,2) = cells_at_HAB.X0(max_soot_cell);
%    max_soot_cells(h,3) = cells_at_HAB.X1(max_soot_cell);
%    
%    max_soot_cells(h,4) = moments(max_soot_cells(h,1)+1,3)/moments(max_soot_cells(h,1)+1,1);
%    max_soot_cells(h,5) = moments(max_soot_cells(h,1)+1,2);
%    max_soot_cells(h,6) = moments(max_soot_cells(h,1)+1,1);
%    
%    ni_tmp = sum(rho(max_soot_cells(h,1)+1)*ni(max_soot_cells(h,1)+1,index_crit:end)'.*dv(index_crit:end));
%    np_tmp = sum(rho(max_soot_cells(h,1)+1)*np(max_soot_cells(h,1)+1,index_crit:end)'.*dv(index_crit:end));
%    
%    nav_tmp = sum(np(max_soot_cells(h,1)+1,:)./ni(max_soot_cells(h,1)+1,:))/m;
%    
% %    nav_agg_only(h) = np_tmp/ni_tmp;
%      nav_agg_only(h) = cells_at_HAB.X1(max_soot_cell)/cells_at_HAB.X0(max_soot_cell);
%    
% %    moments(max_soot_cells(h,1)+1,2);
% %    max_soot_cells(h,4) = moments(cells_at_HAB.CellID,1);
% %    max_soot_cells(h,5) = moments(cells_at_HAB.CellID,2);
%     dp_av = 0;
%     
%     kf = 1.94;
%     Df = 1.8;
% 
% % 
% %     for i = 1:m
% %         np_i = np(max_soot_cell,i);
% %         ni_i = ni(max_soot_cell,i);
% %         nav = np_i/(ni_i + 1*10^-22);
% %         vm = v_m(i);
% %         vp = vm/nav;
% %         
% %         dp_tmp = (6*vp/pi)^(1/3);
% % %         
% % %         dg = dp_tmp*(vm/(vp*kf))^(1/Df);
% % %         
% % %         
% % %         dp = dg*(nav/kf)^Df;
% % %         
% % %         
% % %         dp_av = dp_av + (dp*rho(max_soot_cell)*ni(max_soot_cell,i)*dv(i))...
% % %             /moments(max_soot_cell,1);
% %         
% %     end 
% 
% 
% %--------------------------------------------------------------------------
%     % volume fraction based calculation vp_average = M1/Mp,o
% %--------------------------------------------------------------------------
% if(dp_select == 1)
%     vp = moments(max_soot_cells(h,1)+1,2)/moments(max_soot_cells(h,1)+1,3);
%     max_soot_cells(h,7) = (6*vp/pi)^(1/3);
%     
% %--------------------------------------------------------------------------
%     % surface area based calculation ap_average = SA/Mp,o
% %--------------------------------------------------------------------------
% elseif(dp_select ==2)
%     ap = (cells_at_HAB.Surface_area(max_soot_cell))/moments(max_soot_cells(h,1)+1,3);
%     max_soot_cells(h,7) = (ap/pi)^0.5;
%     
% %--------------------------------------------------------------------------
%     % Volume Fraction to Surface Area  
% %--------------------------------------------------------------------------
% elseif(dp_select == 3)
%     Vn = (moments(max_soot_cells(h,1)+1,2)/moments(max_soot_cells(h,1)+1,1));
%     An = cells_at_HAB.Surface_area(max_soot_cell)/moments(max_soot_cells(h,1)+1,1);
%     max_soot_cells(h,7) = 6*Vn/An;
% 
% 
% 
%     
% end    
%     
% %     max_soot_cells(h,7) = moments(max_soot_cell,1)*moments(max_soot_cell,2)/moments(max_soot_cell,3);
%    
% end

%% plot max annulus line
dp_av_exp = dlmread('./Sun_results/max_dp');% NS flame

% dp_av_exp = dlmread('S_flame_dp_av_MegaDobbin1988_TEM.txt');
plot(dp_av_exp(:,1)/100, dp_av_exp(:,2)/10^9,'ok');

hold on

tol = 0.1;
sp = spap2(20,5,HABs, max_soot_cells(:,7));
fnplt(sp)

axis([0,0.08,0,8*10^-8])
%%
Nav = dlmread('./Sun_results/max_avNpi');

%%
semilogy(Nav(:,1)/100, Nav(:,2),'ok');

%% Nav
hold on
tol = 0.1;
% sp = spaps(HABs, max_soot_cells(:,4),tol*(max(max_soot_cells(:,4)) - min(max_soot_cells(:,4)))/std(max_soot_cells(:,4))   );
% k = 5;
% sp = spapi(optknt(HABs',k),HABs,max_soot_cells(:,4));
% fnplt(sp,2,'g');

sp = spap2(20,5,HABs, max_soot_cells(:,4));
fnplt(sp)



% semilogy(HABs, max_soot_cells(:,4),'r');
% fnplt(sp,2)
axis([0.01,0.06,1,1000])
%%

hold on
tol = 0.1;
% sp = spaps(HABs, max_soot_cells(:,4),tol*(max(max_soot_cells(:,4)) - min(max_soot_cells(:,4)))/std(max_soot_cells(:,4))   );
% k = 5;
% sp = spapi(optknt(HABs',k),HABs,max_soot_cells(:,4));
% fnplt(sp,2,'g');
nav_agg_only(:) = 1.94*(((nav_agg_only(:)*(6/pi)).^(1/3))./max_soot_cells(:,7)).^1.8;

sp = spap2(20,5,HABs, nav_agg_only(:));
fnplt(sp)

% semilogy(HABs, max_soot_cells(:,4),'r');
% fnplt(sp,2)
axis([0.01,0.06,1,1000])

%%
hold on
tol = 0.1;
% sp = spaps(HABs, max_soot_cells(:,4),tol*(max(max_soot_cells(:,4)) - min(max_soot_cells(:,4)))/std(max_soot_cells(:,4))   );
% k = 5;
% sp = spapi(optknt(HABs',k),HABs,max_soot_cells(:,4));
% fnplt(sp,2,'g');

sp = spap2(20,5,HABs, nav_agg_only(:));
fnplt(sp)

% semilogy(HABs, max_soot_cells(:,4),'r');
% fnplt(sp,2)
axis([0.01,0.06,1,1000])

%% centreline

T_centreline = dlmread('./santoro_centreline_exp_results/central_line_T_exp');

plot(T_centreline(:,1)*10^-2,T_centreline(:,2),'o');
%%
hold on 
 
plot(cell_data.cellcentre_0(1:200),cell_data.T(1:200)); 
axis([0,0.1,200,1800])

%% centreline fv

fv_centreline = dlmread('./santoro_centreline_exp_results/central_line_fv_exp');
semilogy(fv_centreline(:,1)*10^-2,fv_centreline(:,2),'o');

%% 
hold on

semilogy(cell_data.cellcentre_0(1:200),moments(1:200,2)*10^6);
axis([0.01,0.08,10^-3,10^2]);
% axis([0,0.1,200,1800])

%% Radial Temp
 
chosen_HAB = [0.003,0.005,0.007,0.01, 0.02, 0.05, 0.069];
all_habs = unique(cell_data.cellcentre_0);

H = size(chosen_HAB,2);

T_radial = zeros(100,2,H);
OH_radial= zeros(100,2,H);
rho_radial= zeros(100,2,H);
acet_radial = zeros(100,2,H);

for h = 1:size(chosen_HAB,2)
   HAB = chosen_HAB(h);
   j=1;
   
   while(all_habs(j)<HAB)
      j=j+1; 
   end

   
   row_beneath = cell_data.cellcentre_0 == all_habs(j-1);
   row_above   = cell_data.cellcentre_0 == all_habs(j);

   
   cells_above = cell_data(row_above,:);
   cells_beneath = cell_data(row_beneath,:);
   
   
   radial_locs = unique(cells_above.cellcentre_1);

   
   for k = 1:size(cells_above)

       
       T_radial(k,1,h)    = radial_locs(k);
       T_radial(k,2,h)    = 0.5*(cells_beneath.T(k) + cells_above.T(k));

       
       OH_radial(k,1,h)   = radial_locs(k);
       OH_radial(k,2,h)   = 0.5*(cells_beneath.OH(k) + cells_above.OH(k));

       
       rho_radial(k,1,h)  = radial_locs(k);
       rho_radial(k,2,h)  = 0.5*(cells_beneath.drho(k) + cells_above.drho(k));

       
       acet_radial(k,1,h) = radial_locs(k);
       acet_radial(k,2,h) = 0.5*(cells_beneath.C2H2(k) + cells_above.C2H2(k));  
   end  
end

%%
T_3mm_exp = dlmread('./Sun_results/T_3mm.csv');
T_5mm_exp    = dlmread('./Sun_results/T_5mm.csv');
T_10mm_exp   = dlmread('./Sun_results/T_10mm.csv');

 
T_20mm_exp = dlmread('./Sun_results/T_20mm.csv');
T_50mm_exp = dlmread('./Sun_results/T_50mm.csv');
T_70mm_exp = dlmread('./Sun_results/T_70mm.csv');

 
OH_7mm_exp = dlmread('./Sun_results/OH_7mmHAB');
OH_70mm_exp = dlmread('./Sun_results/OH_70mmHAB');

 
acet_7mm_exp = dlmread('./Sun_results/C2H2_7mmHAB');
acet_20mm_exp= dlmread('./Sun_results/C2H2_20mmHAB');


%% early T radial profiles
plot(T_3mm_exp(:,1)*10^-3,T_3mm_exp(:,2),'or');
hold on
plot(T_5mm_exp(:,1)*10^-3,T_5mm_exp(:,2),'ob');
plot(T_10mm_exp(:,1)*10^-3,T_10mm_exp(:,2),'ok');

 
hold on
plot(T_radial(:,1,1),T_radial(:,2,1),'-r');
plot(T_radial(:,1,2),T_radial(:,2,2),'-b');
plot(T_radial(:,1,4),T_radial(:,2,4),'k-');

 
%% Later T radial
plot(T_20mm_exp(:,1)*10^-3,T_20mm_exp(:,2),'or');
hold on
plot(T_50mm_exp(:,1)*10^-3,T_50mm_exp(:,2),'ob');
plot(T_70mm_exp(:,1)*10^-3,T_70mm_exp(:,2),'ok');

 
hold on
plot(T_radial(:,1,5),T_radial(:,2,5),'-r');
plot(T_radial(:,1,6),T_radial(:,2,6),'-b');
plot(T_radial(:,1,7),T_radial(:,2,7),'k-');

 
%% OH radial profiles

 
plot(OH_7mm_exp(:,1)*10^-2,OH_7mm_exp(:,2),'or');
hold on
plot(OH_70mm_exp(:,1)*10^-2,OH_70mm_exp(:,2),'ob');

 
hold on
plot(OH_radial(:,1,3),OH_radial(:,2,3),'-r');
plot(OH_radial(:,1,7),OH_radial(:,2,7),'-b');

 
%% acet radial profiles
 
plot(acet_7mm_exp(:,1)*10^-2,acet_7mm_exp(:,2),'ro');
hold on 
plot(acet_20mm_exp(:,1)*10^-2,acet_20mm_exp(:,2),'bo');
 
plot(acet_radial(:,1,3),acet_radial(:,2,3),'r-');
plot(acet_radial(:,1,5),acet_radial(:,2,5),'b-');

%% save matlab objs

if not(isfolder(test_run_name))
    mkdir(test_run_name)
end

save(strcat(strcat('./',test_run_name),'/',test_run_name,n_steps));

%% get psd at cell id

mode = 9; %1 = psd, 2= X0

for c = 1:size(max_soot_cells,1)
    cell_id = max_soot_cells(c,1);
   n = ni(cell_id,:);
   n_p = np(cell_id,:);

   hab = cell_data.cellcentre_0(cell_id +1);
   x1 = cell_data.X1(cell_id);
   x0 = cell_data.X0(cell_id);
   nu_fv = cell_data.Nuc_fv(cell_id);
   ox = cell_data.Oxidation_O2_fv(cell_id) + cell_data.Oxidation_OH_fv(cell_id);
   surface_area = cell_data.Surface_area(cell_id);
   g_fv = cell_data.Growth_fv(cell_id);
   x01 = cell_data.X01(cell_id);
   T = cell_data.T(cell_id);

   if(mode ==1 )
       loglog(v_m,n);
     text(v_m(30),10^40,strcat('HAB = ',num2str(hab)));
     text(v_m(30),10^30,strcat('X1 = ',num2str(x1)));
        text(v_m(30),10^25,strcat('X0 = ',num2str(x0)));
        axis([v_m(1),v_m(end),10^20,10^50])
       x = input(strcat('curr cellid = ',num2str(cell_id)));
   elseif(mode ==2 )
       hold on
       plot(hab,x1,'ob');
   elseif(mode==3)
       hold on
       plot(hab,nu_fv,'ro')
   elseif(mode==4)
       hold on
       plot(hab,ox,'ro')
   elseif(mode==5)
       hold on
       plot(hab,surface_area,'bo')
   elseif(mode==6)
       hold on
       plot(x01,g_fv,'ro')
   elseif(mode==7)
       hold on
       plot(hab,T,'bo')
   
    elseif(mode ==8)
       hold on
       plot(hab,x0,'or-');
        text(hab,x0,num2str(cell_id));
   elseif(mode==9)
       
        loglog(v_m,n,'r');
        hold on
        loglog(v_m,n_p,'b');
        txt = strcmp('HAB =',num2str(hab));
        
        disp(txt);
        disp(num2str(hab));
        
%         text(v_m(30),10^40,num2str(hab));
        
      text(v_m(30),10^40,strcat('HAB = ',num2str(hab)));
      text(v_m(5),1,'NS flame 5,000 K fusing model');
%      text(v_m(30),10^30,strcat('X1 = ',num2str(x1)));
%         text(v_m(30),10^25,strcat('X0 = ',num2str(x0)));
%         axis([v_m(1),v_m(end),10^20,10^50])

       x = input(strcat('curr cellid = ',num2str(cell_id)));
       close all
   end
end














