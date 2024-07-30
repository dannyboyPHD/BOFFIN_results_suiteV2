function [int_SVF] = calc_ISVF(HABs, cell_data, moments, cell_info)

int_SVF = zeros(length(HABs),1);
int_surface_area = zeros(length(HABs),1);
int_ox_rate =zeros(length(HABs),1);
int_OH = zeros(length(HABs),1);
int_O2 = zeros(length(HABs),1);
int_T = zeros(length(HABs),1);
int_oh_fv = zeros(length(HABs),1);
int_o2_fv = zeros(length(HABs),1);

for h = 1:length(HABs)
   HAB = HABs(h);
   
   rowfilter_at_HAB = cell_data.cellcentre_0 == HAB;
   cells_at_HAB = cell_data(rowfilter_at_HAB,:);
   
   cell_ids = cells_at_HAB.CellID;
   tmp_SVF = 0;
   tmp_SA = 0;
   tmp_oxrate = 0;
   tmp_oh = 0;
   tmp_o2 = 0;
   tmp_oh_fv = 0;
   tmp_o2_fv = 0;
   tmp_T = 0;

   for c = 1:length(cell_ids)
      tmp_SVF = tmp_SVF + 0.5*cell_info(cell_ids(c)+1,26)*moments(cell_ids(c)+1,2);
      
      tmp_SA = tmp_SA   + 0.5*cell_info(cell_ids(c)+1,26)*cell_data.Surface_area(cell_ids(c)+1)/cell_data.X0(cell_ids(c)+1);
      
      tmp_oxrate = tmp_oxrate   + 0.5*cell_info(cell_ids(c)+1,26)*(cell_data.Oxidation_O2_fv(cell_ids(c)+1) + cell_data.Oxidation_OH_fv(cell_ids(c)+1));
      tmp_o2_fv = tmp_o2_fv + 0.5*cell_info(cell_ids(c)+1,26)*(cell_data.Oxidation_O2_fv(cell_ids(c)+1));
      tmp_oh_fv = tmp_oh_fv + 0.5*cell_info(cell_ids(c)+1,26)*(cell_data.Oxidation_OH_fv(cell_ids(c)+1));
      tmp_oh = tmp_oh   + 0.5*cell_info(cell_ids(c)+1,26)*( cell_data.drho(cell_ids(c)+1)*cell_data.OH(cell_ids(c)+1));
      tmp_o2 = tmp_o2   + 0.5*cell_info(cell_ids(c)+1,26)*( cell_data.drho(cell_ids(c)+1)*cell_data.O2(cell_ids(c)+1));
      
      tmp_T = tmp_T   + 0.5*cell_info(cell_ids(c)+1,26)*cell_data.T(cell_ids(c)+1 );
   end

   int_SVF(h) = tmp_SVF;
   int_surface_area(h) = tmp_SA;
   int_ox_rate(h) = tmp_oxrate;
   int_OH(h) = tmp_oh;
   int_O2(h) = tmp_o2;
   
   int_oh_fv(h) = tmp_oh_fv;
   int_o2_fv(h) = tmp_o2_fv;
  
   int_T(h) = tmp_T;
   
end

end

