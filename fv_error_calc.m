function [err] = fv_error_calc(fv_exp,z,fv)


    % T_exp 2D array
    err = 0;
    for i = 1:length(fv_exp(:,1))
       i_HAB = HAB2row_no(fv_exp(i,1),z);
       
       
       err = abs(fv(i_HAB) - fv_exp(i,2))/fv_exp(i,2) + err;
        
    end
    
    err = err/length(fv_exp(:,1));

end

