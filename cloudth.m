close all
clear all

%------------------------------------------------------------------------------
% REMARKS: 
% output file is written every 90 s
% to read data every hour, set dt = 40 (90*40=3600s=1hr)
% to read data for continuous time series, set dt = 1
%------------------------------------------------------------------------------
           
% enter range of time (hour) interested:
ti = 0; tf = 12; 
dt = 1;
cmt = ti*3600/90+1:dt:tf*3600/90+1; %time for cm1out_00000t.nc file 

cth = zeros(1,481);  % cloud-top height matrix, stores cth value for each output file

tic

for i = cmt(1):dt:cmt(end)
    
% filename = '/aos/home/mtang/cm1out_000241.nc';
 filename = sprintf('/gs/project/qrw-161-ad/mtang/dx250_mar07/cm1out_%06d.nc',i);

% read data from output files
qc = ncread(filename,'qc');   
qr = ncread(filename,'qr'); 
qi = ncread(filename,'qi');
qs = ncread(filename,'qs');
qg = ncread(filename,'qg');
z = ncread(filename,'z'); 

qh = (qc + qr + qi + qs + qg);   % in kg/kg
    
nk = size(qh,3);

for k = nk:-1:1
    if logical(max(max(1000.*qh(:,:,k))) >= 0.1)   % convert from kg/kg to g/kg
        cth(i) = z(k);
        break;
    end
end

sprintf('cm1out_%06d.nc read and cth is %.3f km',i,cth(i))

%%WRITE cth AS OUTPUT MATRIX FOR FURTHER ANALYSIS
fname = '/gs/project/qrw-161-ad/mtang/dx250_mar07/matrices';
filename = 'cth_dx250_mar07_1.txt';
dlmwrite(fullfile(fname,filename),cth);


end

toc

