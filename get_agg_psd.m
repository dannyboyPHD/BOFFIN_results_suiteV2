function [ni] = get_agg_psd(cell_id, index_crit,m,psd_col_headings,raw)

% locate cellid row

cell_scalars = raw(cell_id,:);
ni = zeros(1,m);


for i = 2:m+1
%     disp(psd_col_headings{i})
    ni(i-1) = cell_scalars.(convertCharsToStrings(psd_col_headings{i}(:)))(1) ;
    
end


end

