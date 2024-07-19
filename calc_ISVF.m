function [int_SVF] = calc_ISVF(HABs, cell_data, moments, cell_info)

int_SVF = zeros(length(HABs),1);

for h = 1:length(HABs)
   HAB = HABs(h);
   
   rowfilter_at_HAB = cell_data.cellcentre_0 == HAB;
   cells_at_HAB = cell_data(rowfilter_at_HAB,:);
   
   cell_ids = cells_at_HAB.CellID;
   tmp_SVF = 0;

   for c = 1:length(cell_ids)
      tmp_SVF = tmp_SVF + 0.5*cell_info(cell_ids(c)+1,26)*moments(cell_ids(c)+1,2);
   end

   int_SVF(h) = tmp_SVF;
   
end

end

