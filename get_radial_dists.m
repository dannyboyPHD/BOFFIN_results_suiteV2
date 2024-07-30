function [T_radial, OH_radial,rho_radial,acet_radial,total_conc_radial,ox_rate_fv] = get_radial_dists(chosen_HAB,all_habs,cell_data)

H = size(chosen_HAB,2);

T_radial = zeros(100,2,H);
OH_radial= zeros(100,2,H);
rho_radial= zeros(100,2,H);
acet_radial = zeros(100,2,H);

ox_rate_fv = zeros(100,2,H);
total_conc_radial = zeros(100,2,H);

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
       
       total_conc(k,1,h) = radial_locs(k);
       P = 0.5*(cells_beneath.P(k) + cells_above.P(k));
       T = 0.5*(cells_beneath.T(k) + cells_above.T(k));
%        P/(8.314 * T)
       total_conc_radial(k,2,h) = 10^5/(8.314 * T);
       
       
       ox_rate_fv(k,1,h) = radial_locs(k);
       ox_rate_fv(k,2,h) = 0.5*(cells_beneath.Oxidation_O2_fv(k) + cells_beneath.Oxidation_OH_fv(k) ...
           + cells_above.Oxidation_O2_fv(k) + cells_above.Oxidation_OH_fv(k));
   end  
end
end

