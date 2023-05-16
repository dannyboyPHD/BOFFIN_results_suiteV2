%% read csv

cell_data = readtable('test_celldata.csv');
grid = readtable('grid_points.csv');


%% grid size
NX = 200;
NY = 100;
NZ = 1;

total_cells = NX*NY*NZ;

cell_info = zeros(total_cells, 26);
offsets = [1,NX+1, NX+2, (NX+1)*(NY+1),(NX+1)*(NY+1)+1, (NX+1)*(NY+1) + (NX+1), (NX+1)*(NY+1) + (NX+2)];
%% add cellcentre cols to cell_data

cellcentre_0 = zeros(total_cells,1);
cellcentre_1 = zeros(total_cells,1);
cellcentre_2 = zeros(total_cells,1);

cell_data = addvars(cell_data,cellcentre_0);
cell_data = addvars(cell_data,cellcentre_1);
cell_data = addvars(cell_data,cellcentre_2);





%%
% row_no = true cellid +1
for cell_id = 1:total_cells

    anchor = [cell_data.StructuredCoordinates_0(cell_id),cell_data.StructuredCoordinates_1(cell_id),cell_data.StructuredCoordinates_2(cell_id)];
    % find corresponding point data

    p1 = (grid.StructuredCoordinates_0 == anchor(1)) & (grid.StructuredCoordinates_1 == anchor(2)) & (grid.StructuredCoordinates_2 == anchor(3));
    p1 = grid(p1,:);

    anchor_ID = p1.PointID;

    cell_info(cell_id,1:3) = [p1.Points_0,p1.Points_1,p1.Points_2];
    


    for o = 1:length(offsets)
        p = grid.PointID == anchor_ID + offsets(o);
        p = grid(p,:);

        cell_info(cell_id,3*o + 1:3*o + 3) = [p.Points_0, p.Points_1,p.Points_2]; 
    end
    
    % calculate cell centres and add new column to table 
    % TO DO (1)
    
    
    
    cell_info(cell_id,25) = cell_info(cell_id,8) - cell_info(cell_id,2);
    cell_info(cell_id,26) = cell_info(cell_id,8)^2 - cell_info(cell_id,2)^2;
    
    [x,y,z] = get_cellcentre(cell_info(cell_id,1:3),cell_info(cell_id,4:6),cell_info(cell_id,7:9),cell_info(cell_id,13:15));
    cell_data.cellcentre_0(cell_id) = x;
    cell_data.cellcentre_1(cell_id) = y;
    cell_data.cellcentre_2(cell_id) = z;
    
end

%%
% PSD parameters
m = 60;
index_crit = 13;


% psds = zeros()
ni = zeros(total_cells,m);
np = zeros(total_cells,m);
moments = zeros(total_cells, 4);
dps = zeros(total_cells,2);

[v,dv,v_m, psd_headings] = get_pbe_grid_info(cell_data);

for cell_id = 1:total_cells
    ni(cell_id,:) = get_agg_psd(cell_id,index_crit,m,psd_headings,cell_data);
    np(cell_id,:) = get_np(cell_id,index_crit,m,psd_headings,cell_data);
  
    moments(cell_id,1) = sum(  ni(cell_id,:)'.*dv(:)         ); % zeroth moments
    moments(cell_id,2) = sum(  ni(cell_id,:)'.*dv(:).*v_m(:) ); % first moments
    moments(cell_id,3) = sum(  np(cell_id,:)'.*dv(:)         );
    
end


%% calc integrated Soot VF

% get x positions: HAB (m)

HABs = unique(cell_data.cellcentre_0);

int_SVF = zeros(length(HABs),1);

for h = 1:length(HABs)
   HAB = HABs(h);
   
   rowfilter_at_HAB = cell_data.cellcentre_0 == HAB;
   cells_at_HAB = cell_data(rowfilter_at_HAB,:);
   
   cell_ids = cells_at_HAB.CellID;
   tmp_SVF = 0;
%  
   for c = 1:length(cell_ids)
      tmp_SVF = tmp_SVF + 0.5*cell_info(cell_ids(c)+1,26)*moments(cell_ids(c)+1,2);
%       cells_at_HAB.X1(c);
      
%       cells_at_HAB.X1(c);
%     ;
   end
%    cells_at_HAB.X1(c)
   int_SVF(h) = tmp_SVF;
   
end

%%
plot(HABs, 2*pi*int_SVF,'r')

%% find annulus of max soot

[max_soot, maxsoot_cell] = max(moments(:,2)); 

radius_of_max_soot = cell_data.cellcentre_1(maxsoot_cell);
% radius_of_max_soot = 0.0026;

% extract all cells with that radius

HABs = unique(cell_data.cellcentre_0);

max_soot_cells = zeros(length(HABs),5);

for h = 1:length(HABs)
   HAB = HABs(h);
   
   rowfilter_at_HAB = cell_data.cellcentre_0 == HAB;
   
%        & cell_data.cellcentre_1 == radius_of_max_soot;
   
   cells_at_HAB = cell_data(rowfilter_at_HAB,:);
   
   [SVF,max_soot_cell] = max(cells_at_HAB.X1);
   
   max_soot_cells(h,1) = cells_at_HAB.CellID(max_soot_cell);
   max_soot_cells(h,2) = cells_at_HAB.X0(max_soot_cell);
   max_soot_cells(h,3) = cells_at_HAB.X1(max_soot_cell);
   
   max_soot_cells(h,4) = moments(max_soot_cell,3);
   max_soot_cells(h,5) = (cells_at_HAB.X1(max_soot_cell)/...
       (moments(max_soot_cell,3)*cells_at_HAB.X0(max_soot_cell))*(6/pi))^(1/3);
   
%    max_soot_cells(h,4) = moments(cells_at_HAB.CellID,1);
%    max_soot_cells(h,5) = moments(cells_at_HAB.CellID,2);
   
   
end

%% plot max annulus line

semilogy(HABs, max_soot_cells(:,2))
% axis([0,0.08,0,10^-5])

semilogy(HABs, max_soot_cells(:,4));




