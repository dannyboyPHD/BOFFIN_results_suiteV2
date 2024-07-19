function [results] = load_Santoro_exp_data(flame)
results = struct;

% cell_data = addvars(cell_data,cellcentre_0);
if (strcmp(flame,'NS'))
    ISVF = dlmread('./Sun_results/fv_integrated');
    results = setfield(results, 'ISVF',ISVF);
    
    dp_av_exp = dlmread('./Sun_results/max_dp');% NS flame
    results = setfield(results, 'dp_av_exp',dp_av_exp);
    
    Nav = dlmread('./Sun_results/max_avNpi');
    results = setfield(results, 'Nav',Nav);
    
    T_centreline = dlmread('./santoro_centreline_exp_results/central_line_T_exp');
    results = setfield(results, 'T_centreline',T_centreline);
    
    fv_centreline = dlmread('./santoro_centreline_exp_results/central_line_fv_exp');
    results = setfield(results, 'fv_centreline',fv_centreline);
    
    T_3mm_exp = dlmread('./Sun_results/T_3mm.csv');
    results = setfield(results, 'T_3mm_exp',T_3mm_exp);
    
    T_5mm_exp    = dlmread('./Sun_results/T_5mm.csv');
    results = setfield(results, 'T_5mm_exp',T_5mm_exp);
    
    T_10mm_exp   = dlmread('./Sun_results/T_10mm.csv');
    results = setfield(results, 'T_10mm_exp',T_10mm_exp);

    T_20mm_exp = dlmread('./Sun_results/T_20mm.csv');
    results = setfield(results, 'T_20mm_exp',T_20mm_exp);
    
    T_50mm_exp = dlmread('./Sun_results/T_50mm.csv');
    results = setfield(results, 'T_50mm_exp',T_50mm_exp);
    
    T_70mm_exp = dlmread('./Sun_results/T_70mm.csv');
    results = setfield(results, 'T_70mm_exp',T_70mm_exp);

    OH_7mm_exp = dlmread('./Sun_results/OH_7mmHAB');
    results = setfield(results, 'OH_7mm_exp',OH_7mm_exp);
    
    OH_70mm_exp = dlmread('./Sun_results/OH_70mmHAB');
    results = setfield(results, 'OH_70mm_exp',OH_70mm_exp);

    acet_7mm_exp = dlmread('./Sun_results/C2H2_7mmHAB');
    results = setfield(results, 'acet_7mm_exp',acet_7mm_exp);
    
    acet_20mm_exp= dlmread('./Sun_results/C2H2_20mmHAB');
    results = setfield(results, 'acet_20mm_exp',acet_20mm_exp);
end

end
