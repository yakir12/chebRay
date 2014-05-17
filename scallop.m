close all
clear all
clc
chebfunpref('factory');
chebfunpref('splitting','on');
%% Let there be light
app = 286;
L = 500;
r = app/2; % aperture radius
maxbeta = atan(r/L); % angle between optical axis and periphery rays
dombeta = [0,maxbeta]; % domain of beta where light is pointing down at -pi/2
beta = chebfun(@(beta) beta,dombeta);
p = 0*beta+1i*L; % a point source at (0,4)
I = exp((beta-pi/2)*1i);
%% Build the refractive and image planes
domx = [-300,300]; % domain of optical setup 
x = chebfun(@(x) x,domx);
cd = 100./(5e-9*x.^4 + 1) -100; % conrnea distal interface
ld = 98./(3.5e-9*x.^4 + 1) -98; % lens distal interface
lp = 0.0012*x.^2; % lens proximal interface
mr = 0.001*x.^2; % mirror
s = [cd,ld,lp,lp,lp,mr]; % quasimatrix
thick = [0,50,250,40,100,180]; % layer thicknesses
s = s-cumsum(thick);
h = s(r,1); % initial distance between the x-axis and cd at the aperture
s = s-h; % translation of the setup 
ri = [1.34,1.37,1.52,1.36,1.38,1.36]; % refractive indices
refractTrue = [true,true,true,true,true,false]; % is the interface refractive or reflective?
%% Ray tracing
ind = [1:6,5,4]; % the order with which the light will travel through the interfaces
direction = [1,diff(ind)] == -1; % the direction the light is approaching the interface
lind = length(ind);
ps = repmat(p,1,lind+1); % initial starting point plus all interfaces
ds = diff(s); % the surface's slope
N = atan(ds)+pi/2; % the surface-normal angle
for i = 1:lind
    p = chebfun(@(beta) intersection(s(:,ind(i)),p(beta),I(beta),x), dombeta,'vectorize'); % finding the intersection points as a function of beta
    if i < lind % don't refract a second time through the distal retina
        if refractTrue(ind(i)) % refract or reflect
            I = refract(N(:,ind(i))+direction(i)*pi,p,I, ri(ind(i)),ri(ind(i)+1)); % find new direction after refraction   
        else
            I = reflect(N(:,ind(i)),p,I); % find new direction after refraction   
        end
    end
    ps(:,i+1) = p;
end
%% Point spread function
sind = [4,5]; % the two image-plane interfaces
psind = [5,9;6,8]; % the two intersection points with both of the image-plane interfaces
l2 = cell(1,2);
w = cell(1,2);
for i = 1:2 % loop through both interfaces of the image-plane 
    loa = cumsum(abs(diff(x+1i*s(:,sind(i))))); % arc length of the image-plane
    loa = loa-loa(0); % adjust arc length to the optical axis
    w{i} = inv2(loa{0,domx(2)}); % arc length from the optical axis as a function of the 'direct' distance to the optical axis
    P2d = r*sqrt(beta/dombeta(2)); % create distances
    d2beta = inv2(real(ps(:,2))); % convert distances to betas
    l2{i} = repmat(chebfun(0,[0,1]),1,2);
    for j = 1:2   
        l = loa(real(ps(:,psind(i,j)))); % interesction of rays expressed as arc length to the optical axis
        l = abs(l); % due to the symmetry
        d2beta = chebfun(@(x) d2beta(x),minandmax(P2d),'splitting','off'); % remove splitting
        l2{i}(:,j) = l(d2beta(P2d)); % redistribute the PSF
    end
end
%% Plot some rays through the system
doubleit = @(x,y) [[flipud(-reshape(x(2:end),[],1));x(:)], [flipud(reshape(y(2:end),[],1));y(:)]];
nbs = 8;
bs = linspace(dombeta(1),dombeta(2),nbs);
fid = fopen('xyscallop.txt','w');
xl = real(ps(bs,:));
y = imag(ps(bs,:));

for i = 2:2:nbs
    xl(i,:) = fliplr(xl(i,:));
end
for i = 2:2:nbs
    y(i,:) = fliplr(y(i,:));
end
mat1 = doubleit(reshape(xl',[],1),reshape(y',[],1))';
mat = num2cell(mat1);
fprintf(fid,'%f %f\n',mat{:});
fclose(fid);

xs = linspace(domx(1),domx(2),1e2)';
ys = s(xs,:);
for i = 1:size(s,2)
    mat1 = [xs,ys(:,i)]';
    fid = fopen(sprintf('xyscallop%u.txt',i),'w');
    mat = num2cell(mat1);
    fprintf(fid,'%f %f\n',mat{:});
    fclose(fid);
end

xs = linspace(domx(1),domx(2),1e2)';
ys = s(xs,:);
for i = 1:size(s,2)-1
    mat1 = [xs,ys(:,i)];
    mat2 = [flipud(xs),flipud(ys(:,i+1))];
    mat3 = [mat1;mat2]';
    fid = fopen(sprintf('xyscallop%uarea.txt',i),'w');
    mat = num2cell(mat3);
    fprintf(fid,'%f %f\n',mat{:});
    fclose(fid);
end
%% Plot the deviation of rays from the optical axis for a scallops eye
for i = 1:2    
    xl2 = linspace(l2{i}.domain(1),l2{i}.domain(2),1e6); % sample across all probabilities
    [psf,psfx] = hist(reshape(l2{i}(xl2,:),[],1),199); % represnt distribution with a histogram
    psf = psf./w{i}(psfx); % adjusting the psf
    % point spread function must sum to one
    mat = doubleit(psfx,psf);
    matpsf = [mat(:,1),mat(:,2)/trapz(mat(:,1),mat(:,2))];
    
    nbs = 1e2;
    bs = linspace(l2{i}.domain(1),l2{i}.domain(2),nbs);
    fid = fopen(sprintf('icdf2scallop%u.txt',i),'w');
    xy = doubleit(rad2deg(bs(:)*dombeta(2)),l2{i}(bs(:),1));
    xy2 = doubleit(rad2deg(bs(:)*dombeta(2)),l2{i}(bs(:),2));
    
    mat = num2cell([xy,xy2(:,2)]');
    fprintf(fid,'x in out\n');
    fprintf(fid,'%f %f %f\n',mat{:});
    fclose(fid);

    fid = fopen(sprintf('psf2scallop%u.txt',i),'w');
    mat = num2cell(matpsf');
    fprintf(fid,'%f %f\n',mat{:});
    fclose(fid);
end