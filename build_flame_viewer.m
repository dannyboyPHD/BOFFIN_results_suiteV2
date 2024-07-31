
flame_viewer(sim_results,disp_names,flame,X,Y,Z);


function flame_viewer(sim_results, disp_names,flame,X,Y,Z)

fig = uifigure;
gl = uigridlayout(fig,[4,3]);
gl.RowHeight = {'fit','fit','1x','1x'};
gl.ColumnWidth = {'1x','1x','1x'};

%menu for save
menu = uimenu(fig);


menu.Text = 'save plots';
menu.MenuSelectedFcn = @menu_save;





% flame plot 
flame_ax = uiaxes(gl);
flame_ax.Layout.Row = [3,4];
flame_ax.Layout.Column = 1;

s = pcolor(flame_ax,X,Y,Z);
set(s, 'EdgeColor', 'none');
colormap(flame_ax,'jet');
s.FaceColor = 'interp';
xlabel(flame_ax,'r (m)', 'interpreter','tex');
ylabel(flame_ax,'HAB (m)', 'interpreter','tex');
axis(flame_ax,[0,0.02,0,0.16]);
hold(flame_ax,'on');

T = evalin('base','load(sim_results{1})');

path_x = T.cell_data.cellcentre_0(T.max_soot_cells(:,1));
path_y = T.cell_data.cellcentre_1(T.max_soot_cells(:,1));
plot(flame_ax,path_y(2:end), path_x(2:end), 'b','LineWidth',2)

[start_cell,i_max_soot] = eta2cell_id(0.36,flame,T.max_soot_cells,T.cell_data);


plot(flame_ax,T.cell_data.cellcentre_1(start_cell), ...
            T.cell_data.cellcentre_0(start_cell), 'ob',...
            'MarkerSize',10, 'MarkerFaceColor', 'b');


        
% PSD aggreate plot
n_ax = uiaxes(gl);
n_ax.Layout.Row = 3;
n_ax.Layout.Column = 2;
%priamry particles
np_ax = uiaxes(gl);
np_ax.Layout.Row = 4;
np_ax.Layout.Column = 2;

text_label = uilabel(gl);
text_label.Layout.Row = [2,4];
text_label.Layout.Column = 3;


eta_text = sprintf('%.2f',0.36);
slider_text = strcat('$\eta = $',eta_text);


slider_title = uilabel(gl);
slider_title.Layout.Row = 1;
slider_title.Layout.Column = 1;
slider_title.Text = slider_text;
slider_title.FontSize = 16;
slider_title.HorizontalAlignment = 'center';
slider_title.Interpreter = 'latex';




output_summary = {};

for r = [1:length(sim_results)]
    
    
    T2 = evalin('base', sprintf('load(sim_results{%i})',r));
    
    [start_cell,i_max_soot] = eta2cell_id(0.36,flame,T2.max_soot_cells,T2.cell_data);
    
    n = T2.ni(start_cell,:);
    n_p = T2.np(start_cell,:);
    disp(disp_names{r});
     
    loglog(n_ax,T2.v_m,n,'-','DisplayName',disp_names{r});
    xlabel(n_ax,'v (m^3)' )
    ylabel(n_ax,'n(v)       (m^-3 m^-3)')
    hold(n_ax,'on');
    
    loglog(np_ax,T2.v_m,n_p,'-','DisplayName',disp_names{r});
    xlabel(np_ax,'v (m^3)' )
    ylabel(np_ax,'n_p(v)     (m^-3 m^-3)')
    hold(np_ax,'on');
    
    
    
    
    x0 = T2.moments(start_cell,1);
    x1 = T2.moments(start_cell,2);
    Nav = T2.max_soot_cells(i_max_soot,4);
    dp_av = T2.max_soot_cells(i_max_soot,7);

    sim_summary = gen_sim_summary(disp_names{r},x0,x1,Nav,dp_av);
    output_summary{end+1} = sim_summary;
           
end
legend(n_ax);
legend(np_ax)
text_label.Text = output_summary;
text_label.VerticalAlignment = 'top'; 





% eta slider
eta_slider = uislider(gl);
eta_slider.Value = 0.35;
eta_slider.Layout.Row = 2;
eta_slider.Layout.Column = 1;
eta_slider.Limits = [0.1,0.6];
eta_slide.MajorTicks = 0:0.5:0.6;
eta_slider.ValueChangingFcn = @etaChangingFcn;
eta_slider.ValueChangedFcn = @etaChangedFcn;

% intiial plotting


    function etaChangingFcn(~,event)
       eta = round(event.Value,2); 
%        c = eta2cell_id(eta,flame,T.max_soot_cells,T.cell_data);
%        hold(flame_ax,'on');
%        flame_ax.Children(1).XData = T.cell_data.cellcentre_1(c);
%        flame_ax.Children(1).YData = T.cell_data.cellcentre_0(c); 
%        hold(flame_ax,'off');
         eta_text = sprintf('%.2f',eta);
         slider_text = strcat('$\eta = $',eta_text);
         slider_title.Text = slider_text;
         
         
         
         
    end


 function etaChangedFcn(~,event)
      eta = event.Value;
      
       [c,d] = eta2cell_id(eta,flame,T.max_soot_cells,T.cell_data);
       hold(flame_ax,'on');
       flame_ax.Children(1).XData = T.cell_data.cellcentre_1(c);
       flame_ax.Children(1).YData = T.cell_data.cellcentre_0(c); 
       hold(flame_ax,'off');
     
     output_summary = {};
     for t = [1:length(sim_results)]
         S = evalin('base', sprintf('load(sim_results{%i})',t));
         
         ni = S.ni;
         np = S.np;
%          max_soot_cells = S.max_soot_cells;
%          cell_data = S.cell_data;
         v_m = S.v_m;
         
         
         [chosen_cell,i] = eta2cell_id(eta,flame,S.max_soot_cells,S.cell_data);
   
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
       
       
    
       x0 = S.moments(chosen_cell,1);
       x1 = S.moments(chosen_cell,2);
       Nav = S.max_soot_cells(i,4);
       dp_av = S.max_soot_cells(i,7);

       sim_summary = gen_sim_summary(disp_names{length(sim_results)-t+1},x0,x1,Nav,dp_av);
       output_summary{end+1} = sim_summary;
     end  
     text_label.Text = output_summary;
       text_label.VerticalAlignment = 'top'; 
 end


% menu save callback
    function menu_save(~,event)
        %n(v)
        SaveName = 'n(v).fig';
        h = figure;
%         h.Visible = 'off';
        
        for i = [1:length(sim_results)]
            x = n_ax.XAxis.Parent.Children(i).XData;
            y = n_ax.XAxis.Parent.Children(i).YData;
            loglog(x,y,'DisplayName',n_ax.Children(i).DisplayName);
            hold on
%             lgndName1 = n_ax.Legend.String{i};
        end
            legend();
%             saveas(h,SaveName,'jpg')
            savefig(h,SaveName)
            delete(h)
 
        SaveName = 'np(v).fig';
        h = figure;
%         h.Visible = 'off';
        for i = [1:length(sim_results)]
            x = np_ax.XAxis.Parent.Children(i).XData;
            y = np_ax.XAxis.Parent.Children(i).YData;
            loglog(x,y,'DisplayName',np_ax.Children(i).DisplayName);
            hold on
%             lgndName1 = n_ax.Legend.String{i};
        end
            legend();
%             saveas(h,SaveName,'jpg')
            savefig(h,SaveName)
            delete(h)
            
        saveas(fig, strcat(flame,'_flameviewer.fig'))
        saveas(fig, strcat(flame,'_flameviewer.jpg','jpg'))
       

        
    end


   



end






