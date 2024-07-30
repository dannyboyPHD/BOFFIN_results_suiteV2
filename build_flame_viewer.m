
flame_viewer(sim_results,disp_names,flame,X,Y,Z);


function flame_viewer(sim_results, disp_names,flame,X,Y,Z)

fig = uifigure;
gl = uigridlayout(fig,[3,3]);
gl.RowHeight = {50,'1x','1x'};
gl.ColumnWidth = {'1x','1x','1x'};

% flame plot 
flame_ax = uiaxes(gl);
flame_ax.Layout.Row = [2,3];
flame_ax.Layout.Column = 1;

s = pcolor(flame_ax,X,Y,Z);
set(s, 'EdgeColor', 'none');
s.FaceColor = 'interp';
xlabel(flame_ax,'r (m)', 'interpreter','tex');
ylabel(flame_ax,'HAB (m)', 'interpreter','tex');
hold(flame_ax,'on');

T = evalin('base','load(sim_results{1})');

path_x = T.cell_data.cellcentre_0(T.max_soot_cells(:,1));
path_y = T.cell_data.cellcentre_1(T.max_soot_cells(:,1));
plot(flame_ax,path_y(2:end), path_x(2:end), 'r','LineWidth',2)

start_cell = eta2cell_id(0.55,flame,T.max_soot_cells,T.cell_data);


plot(flame_ax,T.cell_data.cellcentre_1(start_cell), ...
            T.cell_data.cellcentre_0(start_cell), 'or',...
            'MarkerSize',10, 'MarkerFaceColor', 'r');


        
% PSD aggreate plot
n_ax = uiaxes(gl);
n_ax.Layout.Row = 2;
n_ax.Layout.Column = 2;
%priamry particles
np_ax = uiaxes(gl);
np_ax.Layout.Row = 3;
np_ax.Layout.Column = 2;

for r = [1:length(sim_results)]
    
    T = evalin('base', sprintf('load(sim_results{%i})',r));
    n = T.ni(start_cell,:);
    n_p = T.np(start_cell,:);
    disp(disp_names{r});
     
    loglog(n_ax,T.v_m,n,'-','DisplayName',disp_names{r});
    xlabel(n_ax,'v (m^3)' )
    ylabel(n_ax,'n(v)       (m^-3 m^-3)')
    hold(n_ax,'on')
    
    
    
    loglog(np_ax,T.v_m,n_p,'-','DisplayName',disp_names{r});
    xlabel(np_ax,'v (m^3)' )
    ylabel(np_ax,'n_p(v)     (m^-3 m^-3)')
    hold(np_ax,'on');
       
end
legend(n_ax);
legend(np_ax)






% eta slider
eta_slider = uislider(gl);
eta_slider.Value = 0.36;
eta_slider.Layout.Row = 1;
eta_slider.Layout.Column = 1;
eta_slider.Limits = [0,0.6];
eta_slider.ValueChangingFcn = @etaChangingFcn;
eta_slider.ValueChangedFcn = @etaChangedFcn;

% intiial plotting


    function etaChangingFcn(~,event)
       eta = event.Value; 
       c = eta2cell_id(eta,flame,T.max_soot_cells,T.cell_data);
       hold(flame_ax,'on');
       flame_ax.Children(1).XData = T.cell_data.cellcentre_1(c);
       flame_ax.Children(1).YData = T.cell_data.cellcentre_0(c); 
       hold(flame_ax,'off');
    end


 function etaChangedFcn(~,event)
     for t = [1:length(sim_results)]
         S = evalin('base', sprintf('load(sim_results{%i})',t));
         
         ni = S.ni;
         np = S.np;
         max_soot_cells = S.max_soot_cells;
         cell_data = S.cell_data;
         v_m = S.v_m;
         eta = event.Value;
         
         chosen_cell = eta2cell_id(eta,flame,max_soot_cells,cell_data);
   
         n = ni(chosen_cell,:);
         n_p = np(chosen_cell,:);
%    hab = cell_data.cellcentre_0(chosen_cell +1);
%    x1 = cell_data.X1(chosen_cell);
%    x0 = cell_data.X0(chosen_cell);
%    nu_fv = cell_data.Nuc_fv(chosen_cell);
%    ox = cell_data.Oxidation_O2_fv(chosen_cell) + cell_data.Oxidation_OH_fv(chosen_cell);
%    surface_area = cell_data.Surface_area(chosen_cell);
%    g_fv = cell_data.Growth_fv(chosen_cell);
%    x01 = cell_data.X01(chosen_cell);
%    T = cell_data.T(chosen_cell);
        t = length(sim_results) - t + 1;
       % n(v)
       n_ax.Children(t).XData = v_m;
       n_ax.Children(t).YData = n;
     

       % np(v)
       np_ax.Children(t).XData = v_m;
       np_ax.Children(t).YData = n_p;
     end  
 end


   



end






