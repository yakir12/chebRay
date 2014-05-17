function I = reflect(Ni,p,I)
    Ni = Ni(real(p)); % the angle of the normal as a function of beta
    alpha1  = Ni-angle(I);
    theta = Ni+alpha1-pi;
    I = exp(1i*theta); 