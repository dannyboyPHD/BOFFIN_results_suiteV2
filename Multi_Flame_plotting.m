%% Mulit-Flame Plotting
flame = 'NS';

processed_results_dir = '/Users/danielosullivan/Desktop/FLAME_RESULTS 2/NS_flame/';

sim_results = [{strcat(processed_results_dir, 'NS_matsu1500_sun_growth2760000.mat')}, ...
               {strcat(processed_results_dir, 'NS_tf_ISAFSAF3100000 1.mat')},...
               {strcat(processed_results_dir, 'NS_tfJun27_3075000 1.mat')},...
               {strcat(processed_results_dir,'NS_inf_fusing2525000.mat')}];
           
           
           
disp_names = [{'Matsukawa T_a = 1500 K'}, {'Furnace Black A T_a = 5000 K'},...
    {'Furnace Black B, T_a = 15000 K'},{'Infinite fusing time'}];

experimental_comparison = true;

% Static Plots
integrated_soot_VF = false;
max_soot_path      = true;
centreline         = false;
radial             = false;

%for radial plots
chosen_HAB = [0.003,0.005,0.007,0.01, 0.02, 0.05, 0.069];

% Dynamic PSD Plotting
 PSD_plots         = true;

%% Load Experimental Results

experimental_results = load_Santoro_exp_data(flame);


%% figure creating
if integrated_soot_VF
    isvf_fig = figure;
end

if max_soot_path
    max_soot_fv = figure;
    max_soot_dp_av = figure;
    max_soot_nav = figure;
end

if centreline
   centreline_fv = figure;
   centreline_T   = figure;
end
    
if radial
    radial_T_lower = figure;
    radial_T_higher= figure;
    radial_OH      = figure;
    radial_acet    = figure;
end

%%
for s = [1:length(sim_results)]
    %Load Sim Results
 
        load(sim_results{s});

    if exist('isvf_fig') && integrated_soot_VF
        figure(isvf_fig)
        hold on
        plot(HABs, int_SVF, 'k')
        hold off
         if s == length(sim_results)
            legend(disp_names(:));
        end
    end

    if exist('max_soot_fv') && max_soot_path
        figure(max_soot_fv)

        hold on
        tol = 0.1;
        sp = spap2(20,5,HABs, max_soot_cells(:,3));
        fnplt(sp)
        hold off
         if s == length(sim_results)
            legend(disp_names(:));
        end
    end

    if exist('max_soot_dp_av') && max_soot_path
        figure(max_soot_dp_av)

        hold on
        tol = 0.1;
        sp = spap2(20,5,HABs, max_soot_cells(:,7));
        fnplt(sp)
        axis([0,0.08,0,8*10^-8])
        hold off
         if s == length(sim_results)
            legend(disp_names(:));
        end
    end

    if exist('max_soot_nav') && max_soot_path
        figure(max_soot_nav)

        hold on
        tol = 0.1;
        sp = spap2(20,5,HABs, max_soot_cells(:,4));
        fnplt(sp)
        axis([0.01,0.06,1,1000]);
        hold off
        if s == length(sim_results)
            legend(disp_names(:));
        end
    end

    % Centreline

    if exist('centreline_T') && centreline
        figure(centreline_T)
        hold on
        plot(cell_data.cellcentre_0(1:200),cell_data.T(1:200)); 
        axis([0,0.1,200,1800])
        hold off
    end

    if exist('centreline_fv') && centreline
        figure(centreline_fv)
        if s ==1
            semilogy(cell_data.cellcentre_0(1:200),moments(1:200,2)*10^6);
            axis([0.01,0.08,10^-3,10^2]);
        else
            hold on 
            semilogy(cell_data.cellcentre_0(1:200),moments(1:200,2)*10^6);
            axis([0.01,0.08,10^-3,10^2]);
        end
        
        hold off
    end

    % Radial 

    % T(HABS)
    if exist('radial_T_lower') && radial
        figure(radial_T_lower')
        hold on
        plot(T_radial(:,1,1),T_radial(:,2,1),'-r');
        plot(T_radial(:,1,2),T_radial(:,2,2),'-b');
        plot(T_radial(:,1,4),T_radial(:,2,4),'k-');
        hold off
    end

    if exist('radial_T_higher') && radial
        figure(radial_T_higher')
        hold on
        plot(T_radial(:,1,5),T_radial(:,2,5),'-r');
        plot(T_radial(:,1,6),T_radial(:,2,6),'-b');
        plot(T_radial(:,1,7),T_radial(:,2,7),'k-');
        hold off
    end


    % OH(HABS)
    if exist('radial_OH') && radial
        figure(radial_H')
        hold on
        plot(OH_radial(:,1,3),OH_radial(:,2,3),'-r');
        plot(OH_radial(:,1,7),OH_radial(:,2,7),'-b');
        hold off
    end

    % acet(HABS)
    if exist('radial_acet') && radial
        figure(radial_acet')
        hold on
        plot(acet_radial(:,1,3),acet_radial(:,2,3),'r-');
        plot(acet_radial(:,1,5),acet_radial(:,2,5),'b-');
        hold off
    end
end

%% Experimental Comp




%% Dynamic PSD plotting
eta = 0.55;

Flame_Inter = scatteredInterpolant(cell_data.cellcentre_1, cell_data.cellcentre_0, cell_data.H2O);

x_plot = linspace(min(cell_data.cellcentre_1),max(cell_data.cellcentre_1),100);
y_plot = linspace(min(cell_data.cellcentre_0),max(cell_data.cellcentre_0),200);

[X,Y] = meshgrid(x_plot, y_plot );

Z = Flame_Inter(X,Y);

% chosen_cell = 175;


dynamic = true;

if dynamic
    
    build_flame_viewer;
else
    figure;

for s = [1:length(sim_results)]
    
   load(sim_results{s});
   
   chosen_cell = eta2cell_id(eta,flame,max_soot_cells,cell_data);
    
   n = ni(chosen_cell,:);

  n_p = np(chosen_cell,:);

   hab = cell_data.cellcentre_0(chosen_cell +1);
   x1 = cell_data.X1(chosen_cell);
   x0 = cell_data.X0(chosen_cell);
   nu_fv = cell_data.Nuc_fv(chosen_cell);
   ox = cell_data.Oxidation_O2_fv(chosen_cell) + cell_data.Oxidation_OH_fv(chosen_cell);
   surface_area = cell_data.Surface_area(chosen_cell);
   g_fv = cell_data.Growth_fv(chosen_cell);
   x01 = cell_data.X01(chosen_cell);
   T = cell_data.T(chosen_cell);
   
   path_x = cell_data.cellcentre_0(max_soot_cells(:,1));
   path_y = cell_data.cellcentre_1(max_soot_cells(:,1));
   
   if s ==1
        

        subplot(2,2,2);
        loglog(v_m,n,'-','DisplayName',disp_names{s});
        xlabel('v (m^3)' )
        ylabel('n(v)       (m^-3 m^-3)')
        
        subplot(2,2,4);
        loglog(v_m,n_p,'-','DisplayName',disp_names{s});
        xlabel('v (m^3)' )
        ylabel('n_p(v)     (m^-3 m^-3)')

        
        subplot(2,2,[1,3] );
%         scatter(cell_data.cellcentre_1, cell_data.cellcentre_0, [],cell_data.T)
        s = pcolor(X,Y,Z);
        set(s, 'EdgeColor', 'none');
        s.FaceColor = 'interp';
        xlabel('r (m)', 'interpreter','tex');
        ylabel('HAB (m)', 'interpreter','tex');
%         view(0,90)
        
        hold on
        plot(cell_data.cellcentre_1(chosen_cell), ...
            cell_data.cellcentre_0(chosen_cell), 'or',...
            'MarkerSize',10, 'MarkerFaceColor', 'r');
        plot(path_y(2:end), path_x(2:end), 'r','LineWidth',2)
        hold off


   else     
        subplot(2,2,2);
        hold on 
        loglog(v_m,n,'-','DisplayName',disp_names{s});
        legend();
        hold off
        
        subplot(2,2,4);
        hold on 
        loglog(v_m,n_p,'-','DisplayName',disp_names{s});
       
        legend();
        hold off
       
        subplot(2,2,[1,3]);
        hold on
        plot(cell_data.cellcentre_1(chosen_cell), ...
            cell_data.cellcentre_0(chosen_cell), 'or',...
            'MarkerSize',8, 'MarkerFaceColor', 'r');
        hold off
   end  
end

end


% for s = [1:length(sim_results)]
%     
%    load(sim_results{s});
%    
%    chosen_cell = eta2cell_id(eta,flame,max_soot_cells,cell_data);
%    
%    n = ni(chosen_cell,:);
%    n_p = np(chosen_cell,:);
% 
%    hab = cell_data.cellcentre_0(chosen_cell +1);
%    x1 = cell_data.X1(chosen_cell);
%    x0 = cell_data.X0(chosen_cell);
%    nu_fv = cell_data.Nuc_fv(chosen_cell);
%    ox = cell_data.Oxidation_O2_fv(chosen_cell) + cell_data.Oxidation_OH_fv(chosen_cell);
%    surface_area = cell_data.Surface_area(chosen_cell);
%    g_fv = cell_data.Growth_fv(chosen_cell);
%    x01 = cell_data.X01(chosen_cell);
%    T = cell_data.T(chosen_cell);
%    
%    
%    
%    if s ==1
%        % flame first time
%        path_x = cell_data.cellcentre_0(max_soot_cells(:,1));
%        path_y = cell_data.cellcentre_1(max_soot_cells(:,1));
%        plot(flame_ax,path_y(2:end), path_x(2:end), 'r','LineWidth',2)
%        plot(flame_ax,cell_data.cellcentre_1(chosen_cell), ...
%             cell_data.cellcentre_0(chosen_cell), 'or',...
%             'MarkerSize',10, 'MarkerFaceColor', 'r');
%         
%        hold(flame_ax,'off');
%        
%        % n(v)
%        loglog(n_ax,v_m,n,'-','DisplayName',disp_names{s});
%        xlabel(n_ax,'v (m^3)' )
%        ylabel(n_ax,'n(v)       (m^-3 m^-3)')
%        hold(n_ax,'on')
%        % np(v)
%        loglog(np_ax,v_m,n_p,'-','DisplayName',disp_names{s});
%        xlabel(np_ax,'v (m^3)' )
%        ylabel(np_ax,'n_p(v)     (m^-3 m^-3)')
%        hold(np_ax,'on');
%    else     
%         loglog(n_ax,v_m,n,'-','DisplayName',disp_names{s});
%         legend(n_ax);
%         
%         loglog(np_ax,v_m,n_p,'-','DisplayName',disp_names{s});
%         legend(np_ax);
%    end  
%    
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 






