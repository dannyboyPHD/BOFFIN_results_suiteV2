%% create dir to store raw csvs (cell_Data and grid)

if not(isfolder('./raw_cell_data/'))
    mkdir('./raw_cell_data/')
end

%% create processed results (.mat files)
if not(isfolder('./processed_results/'))
    mkdir('./processed_results/')
end



