function [v,dv,v_m,psd_colheads] = get_pbe_grid_info(raw)


no_pbe_grid_pts = 0;
no_psd_bins = 0;

pbe_v_colheads = {};
psd_colheads = {};

for s = 1:size(raw,2)
    heading = raw.Properties.VariableNames{s};
    disp(heading);
    
   if(length(heading) <2)% heading doesn't concern PBE 
       
   else
      if (strcmp(heading(1:2),'XD'))
        col_head = raw.Properties.VariableNames{s};
        
        pbe_v_colheads{end+1} = col_head;
        
        pbe_cell_no = extractAfter(col_head,'XD');
        pbe_cell_no = str2num(pbe_cell_no);
       
        if(pbe_cell_no > no_pbe_grid_pts)
           no_pbe_grid_pts = pbe_cell_no;
        end
       
      end
      if (strcmp(heading(1:2),'ND'))
          col_head = raw.Properties.VariableNames{s};
        
            psd_colheads{end+1} = col_head;
        
            psd_cell_no = extractAfter(col_head,'ND');
            psd_cell_no = str2num(psd_cell_no);
       
            if(psd_cell_no > no_psd_bins)
            no_psd_bins = psd_cell_no;
         end
          
      end
           
   end
end



% note extra v(0) at nuc size thus:

no_pbe_grid_pts = no_pbe_grid_pts +1;
% NOTE erroreous col ND0000
no_psd_bins = no_psd_bins;


% extract v_i - labelled XD000-XDNNNN
v = zeros(no_pbe_grid_pts,1);

for i = 1:no_pbe_grid_pts
   v(i) = raw.(convertCharsToStrings(pbe_v_colheads{i}(:)))(1);% fixed pbe grid only first HAB needed
end

no_pbe_sections = no_pbe_grid_pts -1;

% create section midpoints
v_m = zeros(no_pbe_sections,1);
dv = zeros(no_pbe_sections,1);

for i = 1:no_pbe_sections
    v_m(i) = (v(i+1) + v(i))*0.5;
    dv(i) = (v(i+1) - v(i));
    
end


end

