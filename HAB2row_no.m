function [i] = HAB2row_no(req_z,z)
cond = true;
i = 1;
while cond
    if(z(i) < req_z)
        i = i+1;
    else
        cond = false;
    end
end

end

