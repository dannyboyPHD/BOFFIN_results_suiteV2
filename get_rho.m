function [rho] = get_rho(cell_id,raw)
cell_scalars = raw(cell_id,:);
rho = cell_scalars.drho; 





end

