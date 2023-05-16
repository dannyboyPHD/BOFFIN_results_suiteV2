function [np] = get_np(cell_id, index_crit,m,psd_col_headings,raw)

cell_scalars = raw(cell_id,:);
np = zeros(1,m);


for i = 2:index_crit+1
%     disp(psd_col_headings{i})
    np(i-1) = cell_scalars.(convertCharsToStrings(psd_col_headings{i}(:)))(1) ;
    
end


for i = 1:m - index_crit
%     disp(psd_col_headings{i})
    off = m + i;
    np(index_crit+i) = cell_scalars.(convertCharsToStrings(psd_col_headings{off}(:)))(1) ;
end



end

