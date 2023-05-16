function [x,y,z] = get_cellcentre(p1,p2,p3,p5)

    x = 0.5*(p1(1) + p2(1));
    y = 0.5*(p2(2) + p3(2));
    z = 0.5*(p1(3) + p5(3));
    
    

end

