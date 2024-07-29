function [cell_data] = calc_add_cell_centres(total_cells,offsets, cell_data,grid,cell_info)

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
    cell_info(cell_id,25) = cell_info(cell_id,8) - cell_info(cell_id,2);
    cell_info(cell_id,26) = cell_info(cell_id,8)^2 - cell_info(cell_id,2)^2;
    
    [x,y,z] = get_cellcentre(cell_info(cell_id,1:3),cell_info(cell_id,4:6),cell_info(cell_id,7:9),cell_info(cell_id,13:15));
    cell_data.cellcentre_0(cell_id) = x;
    cell_data.cellcentre_1(cell_id) = y;
    cell_data.cellcentre_2(cell_id) = z;
    
end


end

