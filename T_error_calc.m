function [err] = T_error_calc(T_exp,z,T)

    % T_exp 2D array
    err = 0;
    for i = 1:length(T_exp(:,1))
       i_HAB = HAB2row_no(T_exp(i,1),z);
       
       
       err = abs(T(i_HAB) - T_exp(i,2))/T_exp(i,2) + err;
        
    end
    
    err = err/length(T_exp(:,1));

end

