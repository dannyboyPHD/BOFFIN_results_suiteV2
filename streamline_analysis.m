%%
test_run_name = 'sinter';

streamline_data = readtable('/home/danny/Desktop/BASE_CASES_SANTORO/NS/streamline_analysis/titania_200core/stream_0.002.csv');
grid = readtable('grid_NS_silica.csv');
%%
cell_data = readtable('celldata_NS_tit_200core_FL0.6.csv');


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
%% calc timescales

% t_agg 
t_agg = zeros(height(streamline_data),1);
t_agg = streamline_data.X0./streamline_data.Aggregation;

% t_nuc
t_nuc = zeros(height(streamline_data),1);
t_nuc = streamline_data.Nuc_fv
% ./(0.5*(streamline_data.XD0001 +streamline_data.XD0002));

% t_sinter

% t_grow
t_nuc = zeros(height(streamline_data),1);
t_grow = streamline_data.Growth_fv
% t_growth
% t_ox

%% plot
for i = 1:height(streamline_data)
   hold on
   plot(cell_data.cellcentre_0(cell_data.CellID == streamline_data.CellID(i)),streamline_data.Growth_fv(i),'o')
    

end

