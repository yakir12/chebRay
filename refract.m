function I = refract(Ni,p,I,n1,n2)
    Ni = Ni(real(p)); % the angle of the normal as a function of beta
    alpha1  = Ni-angle(I);
    alpha2 = asin(n1*sin(alpha1)/n2);
    theta = (Ni+alpha2)-pi;
    I = exp(1i*theta); 