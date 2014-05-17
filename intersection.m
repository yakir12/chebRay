function p = intersection(si,pb,Ib,x)    
    a = imag(Ib)/real(Ib); % the slope of the ray's function
    f = imag(pb)+a*(x-real(pb)); % the ray's function    
    x = roots(si-f);  % the x value where the interface and ray intersect
    y = si(x); % corresponding y value
    p = x+1i*y;