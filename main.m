close all
clear all
clc
chebfunpref('factory');
chebfunpref('splitting','on')
%% Let there be light
app = 1.8;
L = 4;
r = app/2; % aperture radius
maxbeta = atan(r/L); % angle between optical axis and periphery rays
dombeta = [0,maxbeta]; % domain of beta
beta = chebfun(@(beta) beta,dombeta);
p = 0*beta+1i*L; % a point source at (0,4)
I = exp((beta-pi/2)*1i); %  light is pointing down at -90
%% Build the refractive and image planes
domx = [-1,1]; % domain of optical setup 
x = chebfun(@(x) x,domx);
lens1 = -.5*x.^2+.5;
lens2 = .75*x.^2-.75;
retina = -sqrt(3^2-x.^2);
s = [lens1,lens2,retina]; % quasimatrix
h = s(r,1); % initial distance between the x-axis and lens1 at the aperture
s = s-h; % translation of the setup 
ri = [1,1.45,1.3]; % refractive indices
%% Ray tracing
ps = repmat(p,1,3+1); % initial starting point plus 3 interfaces
ds = diff(s); % the surface's slope
N = atan(ds)+pi/2; % the surface-normal angle
for i = 1:3
    p = chebfun(@(beta) intersection(s(:,i),p(beta),I(beta),x), dombeta,'vectorize'); % finding the intersection points as a function of beta
    if i < 3 % don't refract through the retina
        I = refract(N(:,i),p,I,ri(i),ri(i+1)); % find new direction after refraction   
    end
    ps(:,i+1) = p;
end
%% Point spread function
loa = cumsum(abs(diff(x+1i*s(:,end)))); % arc length of the image-plane
loa = loa-loa(0); % adjust arc length to the optical axis
l = loa(real(ps(:,end))); % interesction of rays expressed as arc length to the optical axis
l = abs(l); % due to the symmetry
%% 3D correction for the aperture
P2d = r*sqrt(beta/dombeta(2)); % create distances
d2beta = inv(real(ps(:,2))); % convert distances to betas
l2 = l(d2beta(P2d)); % redistribute the PSF
%% 3D correction for the image-plane
% % % rr = s(:,end); % full retina curve
% % % rr = rr{1e-1,domx(2)}; % retina curve from the optical axis outwards
% % % rr = rr-min(rr); % all values must be non-negative
% % % xx = inv(rr);
% % % % % rr = rr{0,domx(2)};
% % % % % % xx = (rr); % necessary for revolution around the y-axis
% % % % domy = minandmax(rr)';
% % % % % domy = domy-min(domy);
% % % % yy = chebfun('y',domy); %fix this by inversing the circle to rotatearound the x-axis: x should range from 0 to 0.1716 and y will range from 0 to 1.
% % % % xx = sqrt(3^2-yy.^2);
% % % sr = 2*pi*cumsum(xx.*sqrt(1+abs(diff(xx)).^2)); % surface of revolution of the retina 
% % % ffdgfdgfdgd

w = inv(loa{0,domx(2)}); % arc length from the optical axis as a function of the 'direct' distance to the optical axis
%% Plot 2D point spread function
doubleit = @(x,y) [[flipud(-reshape(x(2:end),[],1));x(:)], [flipud(reshape(y(2:end),[],1));y(:)]]; % help function to double the two sides
xl = linspace(dombeta(1),dombeta(2),1e6); % sample across beta's domain
[psf,psfx] = hist(l(xl),100); % represnt distribution with a histogram
psfxy = doubleit(psfx,psf); % rotational symmetry
psfxy(:,2) = psfxy(:,2)/trapz(psfxy(:,1),psfxy(:,2)); % point spread function must sum to one
%% Plot rays through the system
nbs = 8; % plot 8 rays through the system
bs = linspace(dombeta(1),dombeta(2),nbs); % 8 discrete rays
fid = fopen('xy.txt','w');
x = real(ps(bs,:));
y = imag(ps(bs,:));

for i = 2:2:nbs
    x(i,:) = fliplr(x(i,:));
end
for i = 2:2:nbs
    y(i,:) = fliplr(y(i,:));
end
mat1 = doubleit(reshape(x',[],1),reshape(y',[],1))';
mat = num2cell(mat1);
fprintf(fid,'%f %f\n',mat{:});
fclose(fid);

nbs = 100; 
bs = linspace(l.domain(1),l.domain(2),nbs);
fid = fopen('icdf.txt','w');
mat = num2cell(doubleit(rad2deg(bs),l(bs))');
fprintf(fid,'%f %f\n',mat{:});
fclose(fid);

fid = fopen('psf.txt','w');
mat = num2cell(psfxy');
fprintf(fid,'%f %f\n',mat{:});
fclose(fid);
%% Plot 3D point spread function
xl2 = linspace(l2.domain(1),l2.domain(2),1e6); % sample across all probabilities
[psf,psfx] = hist(l2(xl2),99); % represnt distribution with a histogram
psf = psf./w(psfx); % adjusting the psf
psfxy = doubleit(psfx,psf); % rotational symmetry
psfxy(:,2) = psfxy(:,2)/trapz(psfxy(:,1),psfxy(:,2)); % point spread function must sum to one

nbs = 1e2;
bs = linspace(l2.domain(1),l2.domain(2),nbs);
fid = fopen('icdf2.txt','w');
mat = num2cell(doubleit(rad2deg(bs),l2(bs))');
fprintf(fid,'%f %f\n',mat{:});
fclose(fid);

fid = fopen('psf2.txt','w');
mat = num2cell(psfxy');
fprintf(fid,'%f %f\n',mat{:});
fclose(fid);