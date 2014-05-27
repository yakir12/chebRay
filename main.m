close all
clear all
clc
chebfunpref('factory');
chebfunpref('splitting', 'on')
%% Let there be light
app = 1.8;
L = 4;
r = app/2; % aperture radius
maxbeta = atan(r/L); % angle between optical axis and periphery rays
dombeta = [0, maxbeta]; % domain of beta
beta = chebfun(@(beta) beta, dombeta);
p = 0*beta+1i*L; % a point source at (0, 4)
I = exp((beta-pi/2)*1i); %  light is pointing down at -90
%% Build the refractive and image planes
domx = [-1, 1]; % domain of optical setup 
x = chebfun(@(x) x, domx);
lens1 = -.5*x.^2+.5;
lens2 = .75*x.^2-.75;
retina = -sqrt(3^2-x.^2);
z = [lens1, lens2, retina]; % quasimatrix
h = z(r, 1); % initial distance between the x-axis and lens1 at the aperture
z = z-h; % translation of the setup 
ri = [1, 1.45, 1.3]; % refractive indices
%% Ray tracing
ps = repmat(p, 1, 3+1); % initial starting point plus 3 interfaces
dz = diff(z); % the surface's slope
N = atan(dz)+pi/2; % the surface-normal angle
for i = 1:3
    p = chebfun(@(beta) intersection(z(:, i), p(beta), I(beta), x),  dombeta, 'vectorize'); % finding the intersection points as a function of beta
    if i < 3 % don't refract through the retina
        I = refract(N(:, i), p, I, ri(i), ri(i+1)); % find new direction after refraction   
    end
    ps(:, i+1) = p;
end
%% Point spread function
loa = cumsum(abs(diff(x+1i*z(:, end)))); % arc length of the image-plane
loa = loa-loa(0); % adjust arc length to the optical axis
l = loa(real(ps(:, end))); % interesction of rays expressed as arc length to the optical axis
l = abs(l); % due to the symmetry
%% 3D correction for the aperture
P2d = r*sqrt(beta/dombeta(2)); % create distances
d2beta = inv(real(ps(:, 2))); % convert distances to betas
l2 = l(d2beta(P2d)); % redistribute the PSF
%% 3D correction for the image-plane
w = inv(loa{0, domx(2)}); % arc length from the optical axis as a function of the 'direct' distance to the optical axis
%% Plot rays through the system
nbs = 8; % plot 8 rays through the system
bs = linspace(dombeta(1),dombeta(2),nbs); % 8 discrete rays
x = real(ps(bs,:));
y = imag(ps(bs,:));
x = [-x(2:end,:);x]';
y = [y(2:end,:);y]';

figure
plot(z,'b')
hold on
plot(x,y,'k')
title('A simple case study')
axis equal
%% Plot 2D point spread function
doubleit = @(x,y) [[flipud(-reshape(x(2:end),[],1));x(:)], [flipud(reshape(y(2:end),[],1));y(:)]]; % help function to double the two sides
xl = linspace(dombeta(1),dombeta(2),1e6); % sample across beta's domain
[psf,psfx] = hist(l(xl),100); % represnt distribution with a histogram
psfxy = doubleit(psfx,psf); % rotational symmetry
psfxy(:,2) = psfxy(:,2)/trapz(psfxy(:,1),psfxy(:,2)); % point spread function must sum to one

nbs = 100; 
bs = linspace(l.domain(1),l.domain(2),nbs);
mat = doubleit(rad2deg(bs),l(bs));
figure 
plot(mat(:,1),mat(:,2),'k')
xlabel('\beta (\circ)')
ylabel('Deviation')
title('Deviations of rays from the optical axis as a function of \beta (2D)')
%
figure 
semilogy(psfxy(:,1),psfxy(:,2),'k')
xlabel('Deviation')
ylabel('Probability')
title('Point spread functions of the optical system (2D)')
%% Plot 3D point spread function
xl2 = linspace(l2.domain(1),l2.domain(2),1e6); % sample across all probabilities
[psf,psfx] = hist(l2(xl2),99); % represnt distribution with a histogram
psf = psf./w(psfx); % adjusting the psf
psfxy = doubleit(psfx,psf); % rotational symmetry
psfxy(:,2) = psfxy(:,2)/trapz(psfxy(:,1),psfxy(:,2)); % point spread function must sum to one

nbs = 100; 
bs = linspace(l2.domain(1),l2.domain(2),nbs);
mat = doubleit(rad2deg(bs),l2(bs));
figure 
plot(mat(:,1),mat(:,2),'r')
xlabel('\beta (\circ)')
ylabel('Deviation')
title('Deviations of rays from the optical axis as a function of \beta (3D)')

figure 
semilogy(psfxy(:,1),psfxy(:,2),'r')
xlabel('Deviation')
ylabel('Probability')
title('Point spread functions of the optical system (3D)')