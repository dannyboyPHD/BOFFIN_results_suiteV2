function [cell_id,row_max_soot] = eta2cell_id_by_path(eta,flame, pathline_cellids, celldata)
if(strcmp(flame,'NS'))
    vol_flow_rate = 3.85; % cm^3/s
elseif(strcmp(flame,'IS'))
    vol_flow_rate = 4.60; % cm^3/s
elseif(strcmp(flame,'S'))
    vol_flow_rate = 4.90; % cm^3/s
end

D = 0.156; % cm^2/s diffusion co-ef ethylene-Nitrogen --> Hirschfelder, J. O., Curtiss, C. F., and Bird, R. B.,
%Molecular Theory o f Gases and Liquids, John Wiley,
%New York, 1954, p. 579.  ---> used in original Santoro paper

z = vol_flow_rate*eta/D; % in cm
z = z/100 % cm --> m

habs = unique(celldata.cellcentre_0);

cells_on_pathline = celldata(ismember(celldata.CellID,pathline_cellids(:,1)),:);
j = 1;
% find cell id
while(habs(j)< z)
    j=j+1; 
end

rowfilter = cells_on_pathline.cellcentre_0 == habs(j);
%row in cell_data = cell_id +1
cell_id = cells_on_pathline(rowfilter,:).CellID;

[q,row_max_soot] = ismember(cell_id,pathline_cellids(:,1),'rows');

   
end

