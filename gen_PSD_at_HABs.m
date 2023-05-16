function [v_m,psds,v] = gen_PSD_at_HABs(raw, z,mode)
no_HABs = size(raw,1);

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

if(length(z)~= no_HABs)
   disp('error ni raw_data, different no of rows to z'); 
end


% extract v_i - labelled XD000-XDNNNN
v = zeros(no_pbe_grid_pts,1);

for i = 1:no_pbe_grid_pts
   v(i) = raw.(convertCharsToStrings(pbe_v_colheads{i}(:)))(1);% fixed pbe grid only first HAB needed
end

no_pbe_sections = no_pbe_grid_pts -1;

% create section midpoints
v_m = zeros(no_pbe_sections,1);


for i = 1:no_pbe_sections
    v_m(i) = (v(i+1) + v(i))*0.5;
end

% mode for 1 PBE and 2PBE cases i.e. aggregates and primary particle dists

if(strcmp(mode,'2PBE')) % need to add 1PBE mode
    psd_primary = zeros(no_HABs, no_pbe_sections);
    % no of aggregate bins above cut-off
    i_crit = -no_psd_bins + 2*no_pbe_sections;
    
    psd_agg = zeros(no_HABs,no_pbe_sections - i_crit);
    
    
% skip ND0000 
    for h = 1:no_HABs
        disp(h);
        for i = 1:no_pbe_sections
            psd_primary(h,i) = raw.(psd_colheads{i+1})(h);
        end
        
        for j = no_pbe_sections+1:no_psd_bins
            i_off = j - no_pbe_sections;
            psd_agg(h,i_off) = raw.(psd_colheads{j+1})(h);
        end
        
    end
    
    
    psds = zeros(no_HABs,no_pbe_sections,2);
    psds(:,:,1) = psd_primary;
    psds(:,i_crit+1:no_pbe_sections,2) = psd_agg;
end





end

