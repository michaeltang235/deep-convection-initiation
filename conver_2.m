close all;
clear all;

tic

% enter range of time (hour) interested:
ti = 0; tf = 12; 
dt = 1;
cmt = ti*3600/90+1:dt:tf*3600/90+1; %time for cm1out_00000t.nc file 

% enter range on x-axis 
xmin = -0.5;  % in km
xmax = 0.5;

% initialize average surface convergence, convm, matrix
convm = zeros(1,481); 

for i = cmt(1):dt:cmt(end)

% filename = '/aos/home/mtang/cm1out_000241.nc';
 filename = sprintf('/gs/project/qrw-161-ad/mtang/dx250_mar07/cm1out_%06d.nc',i);

% read data from output files

% xh = ncread(filename,'xh');    % xh in km
xf = ncread(filename,'xf');    % xf in km
% yh = ncread(filename,'yh');    % yh in km
yf = ncread(filename,'yf');    % yf in km
u = ncread(filename,'u');      % u in m/s

nj = size(u,2);

xii = find(xf == xmin);  % index of xmin in xf
xfi = find(xf == xmax);   % index of xmax in xf

delx = 1000*(xf(xfi) - xf(xii));  % convert from km to m
dely = 1000*(yf(end) - yf(1));      

uxfi = 0;
uxii = 0;

dyf = 1000*(yf(2:end) - yf(1:end-1));  % convert from km to m

for j = 1:nj
    uxfi = uxfi + u(xfi,j,1)*dyf(j);  % compute sum of u at xf
    uxii = uxii + u(xii,j,1)*dyf(j);   % compute sum of u at xi
end

convm(i) = (-1/(dely*delx))*(uxfi - uxii);

sprintf('cm1out_%06d.nc read and convm is %.4f x10^-3',i,convm(i)*10^3)
%sprintf('cm1out_%06d.nc read and convm is %.4f',i,convm(i))

%%WRITE  AS OUTPUT MATRIX FOR FURTHER ANALYSIS
fname = '/gs/project/qrw-161-ad/mtang/dx250_mar07/matrices';
filename = 'convm_dx250_mar07_2.txt';
dlmwrite(fullfile(fname,filename),convm);

end

toc
