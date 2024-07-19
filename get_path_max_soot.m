function [max_soot_cells,nav_agg_only] = get_path_max_soot(index_crit,HABs,m,ni,np,rho, dv, cell_data,moments,dp_select)

max_soot_cells = zeros(length(HABs),7);

nav_agg_only = zeros(length(HABs),1);

for h = 1:length(HABs)
   HAB = HABs(h);
   
   rowfilter_at_HAB = cell_data.cellcentre_0 == HAB;
   
   cells_at_HAB = cell_data(rowfilter_at_HAB,:);
   
   [SVF,max_soot_cell] = max(cells_at_HAB.X1);
   
   max_soot_cells(h,1) = cells_at_HAB.CellID(max_soot_cell);
   max_soot_cells(h,2) = cells_at_HAB.X0(max_soot_cell);
   max_soot_cells(h,3) = cells_at_HAB.X1(max_soot_cell);
   
   max_soot_cells(h,4) = moments(max_soot_cells(h,1)+1,3)/moments(max_soot_cells(h,1)+1,1);
   max_soot_cells(h,5) = moments(max_soot_cells(h,1)+1,2);
   max_soot_cells(h,6) = moments(max_soot_cells(h,1)+1,1);
   
   ni_tmp = sum(rho(max_soot_cells(h,1)+1)*ni(max_soot_cells(h,1)+1,index_crit:end)'.*dv(index_crit:end));
   np_tmp = sum(rho(max_soot_cells(h,1)+1)*np(max_soot_cells(h,1)+1,index_crit:end)'.*dv(index_crit:end));
   
   nav_tmp = sum(np(max_soot_cells(h,1)+1,:)./ni(max_soot_cells(h,1)+1,:))/m;
   
   nav_agg_only(h) = np_tmp/ni_tmp;
%      nav_agg_only(h) = cells_at_HAB.X1(max_soot_cell)/cells_at_HAB.X0(max_soot_cell);
   
    dp_av = 0;
    % aggregate morphology assumptions
    kf = 1.94;
    Df = 1.8;

    %--------------------------------------------------------------------------
    % volume fraction based calculation vp_average = M1/Mp,o
    %--------------------------------------------------------------------------
    if(dp_select == 1)
        vp = moments(max_soot_cells(h,1)+1,2)/moments(max_soot_cells(h,1)+1,3);
        max_soot_cells(h,7) = (6*vp/pi)^(1/3);
    %--------------------------------------------------------------------------
        % surface area based calculation ap_average = SA/Mp,o
    %--------------------------------------------------------------------------
    elseif(dp_select ==2)
        ap = (cells_at_HAB.Surface_area(max_soot_cell))/moments(max_soot_cells(h,1)+1,3);
        max_soot_cells(h,7) = (ap/pi)^0.5;
    %--------------------------------------------------------------------------
        % Volume Fraction to Surface Area  
    %--------------------------------------------------------------------------
    elseif(dp_select == 3)
        Vn = (moments(max_soot_cells(h,1)+1,2)/moments(max_soot_cells(h,1)+1,1));
        An = cells_at_HAB.Surface_area(max_soot_cell)/moments(max_soot_cells(h,1)+1,1);
        max_soot_cells(h,7) = 6*Vn/An;    
    end    
end




end

