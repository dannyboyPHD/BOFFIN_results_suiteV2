function [dp_i] = calc_dp_i(nav_i,vm_i)

vp = vm_i/nav_i;
dp_i = (6*vp/pi)^(1/3);

end

