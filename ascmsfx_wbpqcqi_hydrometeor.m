close all
clear all

%%----------------------------------------------------------------------------
% REMARKS:
% If cm1out files are written every 90s, set dcmt=90;
% If cm1out files are written evbery 180s, set dcmt=180;
% This code generates profile of ascending mass flux every hour
% and calculates lcl in the range of x-location interested
% This code generates amf arrays which store the following quantities in
% the chosen region
% (i) total ascending mass flux, (ii) total area of ascending points
% (iii) avg. vertical velocity of asc. points
% (iv) avg. buoyancy of asc. points, (v) avg. qh (hydro meteor.) of asc. points, (vi) mean LCL over the hour
% Constraints for asc. points: wh>0 (w-positive), b>0 (buoyancy-positive),
% and qc+qi>0.1 g/kg

% note: this code uses initial state theta_rho to calculate buoyancy, which
% is different from other versions
%%-----------------------------------------------------------------------------
tic

% initialize array for storing y-avg. cloud mass flux matrix for each hour
asmf = { };

for seq = 6:6  %((ti,tf) = (0,1)) 1:12 for all hours
    
% enter range of time (hour) interested:
% ti = 0; tf = 12;
ti = 0+(seq-1)*1; tf = 1+(seq-1)*1;

dt = 4; % time step for reading cm1 output files (make sure number of files read is the same for all cases)
dcmt = 90; % time step between each cm1out file written (s)

% set range of x intersted:
xmin = -60; % initial x-pos. in km
xmax = 60; % final x-pos. in km

% directory to read output files (cm1 files) from:
fname = '/home/tang235/projects/ctb-kirshbau/tang235/dx250_mar07';   % beluga
% fname = '/gs/project/qrw-161-ad/mtang/dx250_mar07';   % guillimin
% fname = '/storage/mtang/dx250_mar07';   % windmill

% % directory for storing msfx matrix and name of matrix:
fname1 = '/home/tang235/projects/ctb-kirshbau/tang235/dx250_mar07/matrices';
% fname1 = '/gs/project/qrw-161-ad/mtang/dx250_mar07/matrices';
% fname1 = '/storage/mtang/dx250_mar07/matrices';
filename1 = sprintf('ascmf_dx250_mar07_wbpqcqi_hydrometeor.mat');   % ascending mass flux matrix

% enter 1 for saving output matrices and 0 for not saving outputs
saveop = 1;

%%---------------------------------------------------------------------------
%%---------------------------------------------------------------------------

%----------------------------------------------------------------------------
% Calculate initial state density potential temperature (for buoyancy)

% INPUT CONSTANTS FROM constants.F
g = 9.81; % m/s^2
Rd = 287.04;
Rv = 461.5;
epsi = Rd/Rv;

% find r_o basic state density from initial file
filename = sprintf('cm1out_%06d.nc',1);
cmfile = fullfile(fname,filename);

qv = ncread(cmfile,'qv');    % qv, water vapor in kg/kg
qc = ncread(cmfile,'qc');    % qc, cloud water in kg/kg
qr = ncread(cmfile,'qr');    % qr, rain 
qs = ncread(cmfile,'qs');    % qs, snow
qi = ncread(cmfile,'qi');    % qi, cloud ice in kg/kg
qg = ncread(cmfile,'qg');    % qg, graupel 

qt = qv + qc + qr + qs + qi + qg; % total mixing ratio

th = ncread(cmfile,'th'); % potential temp. in K

% get half-level size
ni=size(th,1);
nj=size(th,2);
nk=size(th,3);

% calculate density potential temperature
thr = th.*(1+qv./epsi)./(1+qt);  % theta_r

% calculate horizontal mean of initial state density potential temperature
% after reshaping each layer of thr into column vector
thr0 = zeros(1,1,nk);
for k = 1:nk
    thr0(k) = mean(reshape(thr(:,:,k),[ni*nj 1]));  
end

% repmat thr0 into 3-d matrix
thr03d = repmat(thr0, [ni nj 1]);   % initial state hori. avg. (whole domain) density potential temp. vs z

% END Calculating initial state density potential temperature
%----------------------------------------------------------------------------

% field name for structure 
fieldname = sprintf('f%d%d',ti,tf);

cmt = ti*3600/dcmt+1:dt:tf*3600/dcmt+1; %time for cm1out_00000t.nc file 

% INPUT CONSTANTS FROM constants.F
g = 9.81; % m/s^2
p00 = 1.0e5;
Rd = 287.04;
Cp = 1005.7;
Rv = 461.5;
epsi = Rd/Rv;
Lv = 2501000.0;

for s = cmt(1):dt:cmt(end)

filename = sprintf('cm1out_%06d.nc',s);
cmfile = fullfile(fname,filename);

% read data from output files
xf = ncread(cmfile,'xf');   % xf in km
xh = ncread(cmfile,'xh');   % xh in km
yf = ncread(cmfile,'yf');   % yf in km
z = ncread(cmfile,'z');   % z in km
w = ncread(cmfile,'w');      % w in m/s

qv = ncread(cmfile,'qv');    % qv, water vapor in kg/kg
qc = ncread(cmfile,'qc');    % qc, cloud water in kg/kg
qr = ncread(cmfile,'qr');    % qr, rain 
qs = ncread(cmfile,'qs');    % qs, snow
qi = ncread(cmfile,'qi');    % qi, cloud ice in kg/kg
qg = ncread(cmfile,'qg');    % qg, graupel 

qt = qv + qc + qr + qs + qi + qg; % total mixing ratio
qh = qt - qv;   % total hydrometeor mixing ratio

prs = ncread(cmfile,'prs'); % pressure in Pa
th = ncread(cmfile,'th'); % potential temp. in K

% get half-level size
ni=size(qc,1);
nj=size(qc,2);
nk=size(qc,3);

% get sum of qc and qi
qcqi = qc + qi;

% set threshold qth
qth = 0.1;  % threshold qc+qi in g/kg  
qth = qth/1000;  % qth in kg/kg

% calculate density using density temperature
exn = (prs./p00).^(Rd/Cp);       %exner function
temp = exn.*th;
tempr = temp.*(1+qv./epsi)./(1+qt);     %density temp.

rho = prs./(Rd*(tempr));       %density

% calculate buoyancy using theta_r
thr = th.*(1+qv./epsi)./(1+qt);  % theta_r
% thrm = mean(mean(thr,1),2);  % mean theta_r
% thrm3d = repmat(thrm,[ni nj 1]);

b = g*((thr - thr03d)./thr03d);    % buoyancy (relative to initial state theta_rho) 
% b = g*((thr - thrm3d)./thrm3d);    % buoyancy 

% GET w ONTO HALF LEVEL (z)
nkw = size(w,3);
wh = 0.5*(w(:,:,1:nkw-1) + w(:,:,2:nkw));

% get dxh
dxh = xh(2:end) - xh(1:end-1); % in km

% get indices of xhi and xhf on xhf
xhii = max(find(abs(xh-xmin) < dxh(1)/1.5));
xhfi = min(find(abs(xh-xmax) < dxh(1)/1.5));

% if min is used for xhii and max is used for xhfi, total area calcuated
% would show different dependence on grid spacings


% FIND LCL AT X-LOCATION INTERESTED

% convert z from km to m for parcel_free.m
zpf = 1000.*z;  % in m

% get mean (y-avg.) sounding at each x-location
% calculate y-avg. pressure, temp, and qv along x-axis
prsm = mean(prs,2);
tempm = mean(temp,2);
qvm = mean(qv,2);

% reshape prsm, tempp, and qvm into 2-d matrices as inputs
% for parcel_free.m
p_in = reshape(prsm,[ni nk]);
t_in = reshape(tempm,[ni nk]);
qv_in = reshape(qvm,[ni nk]);

% initialize LCL matrix at the beginning of time
if s == cmt(1)
    lclx = zeros(1,xhfi-xhii+1); % in m
end

% get LCL at each point on x-axis using parcel_free.m
% lclx = zeros(1,xhfi-xhii+1); % in m
for i = xhii:xhfi
    [~,~,lcl,~,~,~,~]=parcel_free(zpf,p_in(i,:),t_in(i,:),qv_in(i,:));
    lclx(i-xhii+1) = lclx(i-xhii+1) + lcl;
sprintf('cm1out_000%d.nc read, ti = %d, tf = %d \n lclx at x = %d found',s,ti,tf,xh(i));
end

% get average LCL at the end of time
if s == cmt(end)
lclxm = zeros(1,xhfi-xhii+1); % in m
lclxm = lclx./(length(cmt));
end


% % find index of lcl at each x in z matrix
% zid = zeros(1,length(lclx));
% for i = 1:length(lclx)
%     zid(i) = find(abs(z - lclx(i)/1000) < 0.09);   % convert lclx from m to km
% end


% CALCULATE ASCENDING MASS FLUX AT EACH X AND Z

% get dxf and dyf
dxf = 1000.*(xf(2:end) - xf(1:end-1));   % in m
dyf = 1000.*(yf(2:end) - yf(1:end-1));   % in m

% find indices of xmin and xmax on xf and yf
xfii = find(xf == xmin);
xffi = find(xf == xmax);

% initialize ascending mass flux matrix at the beginning of time
if s == cmt(1)
    mf = zeros(1,nk);   % total ascending mass flux at each z
    fc = zeros(xhfi-xhii+1,nj,nk);  % freq. counter for ascending points in the chosen region
    aa = zeros(xhfi-xhii+1,nj,nk);  % area of each ascending point in the chosen region
    wa = zeros(xhfi-xhii+1,nj,nk);  % vertical velocity of each acending point in the chosen region
    ba = zeros(xhfi-xhii+1,nj,nk);  % buoyancy of ascending point in the chosen region
    qha = zeros(xhfi-xhii+1,nj,nk);  % qh (hydrometeor) of ascending point in the chosen region
end

%initialize 3-d mass flux matrix which stores mf value at each point
dmf = zeros(xhfi-xhii+1,nj,nk);   % ascending mass flux at each point 

for i = xhii:xhfi
    for j = 1:nj
        for k = 1:nk
            if wh(i,j,k) > 0 && b(i,j,k) > 0 && qcqi(i,j,k) > qth % condition to be ascending point (wh>0, b>0, qc+qi>0.1 g/kg)
                dmf(i-xhii+1,j,k) = rho(i,j,k)*wh(i,j,k)*dxf(i)*dyf(j);
                fc(i-xhii+1,j,k) = fc(i-xhii+1,j,k) + 1;  % incrmenet counter by 1 for the point considered as ascending point
                aa(i-xhii+1,j,k) = fc(i-xhii+1,j,k)*dxf(i)*dyf(j);
                wa(i-xhii+1,j,k) = wa(i-xhii+1,j,k) + wh(i,j,k);
                ba(i-xhii+1,j,k) = ba(i-xhii+1,j,k) + b(i,j,k);
                qha(i-xhii+1,j,k) = qha(i-xhii+1,j,k) + qh(i,j,k);
            end
        end
    end
end

% get total mass flux of all ascending points at each level
for k = 1:nk
    mf(1,k) = mf(1,k) + sum(sum(dmf(:,:,k)));
end

sprintf('cm1out_000%d.nc read, ti = %d, tf = %d \n max. mf = %d, max. fc = %d',s,ti,tf,max(mf),max(max(max((fc)))))

end % end s = cmt


% calculate profiles of n_a, a_a, w_am, b_am, and qh_am:
na = zeros(1,nk);    % total number of ascending points at each z-level
aat = zeros(1,nk);   % total area of ascendng points at each level
wam = zeros(1,nk);   % avg. vertical velocity of ascending points at each level
bam = zeros(1,nk);  % profile of avg. b_c vs z
qham = zeros(1,nk);   % avg. qh of ascending points at each level

% get total number of ascending points at each z-level
for k = 1:nk
    na(k) = sum(sum(fc(:,:,k),1),2);
end

% get total area of ascending points at each level
for k = 1:nk
    aat(k) = sum(sum(aa(:,:,k),1),2);
end

% calculate avg. vertical velocity of ascending points at each z level
for k = 1:nk
    if na(k) > 0
        wam(k) = sum(sum(wa(:,:,k),1),2)/na(k);
    end
end

% calculate avg. buoyancy of ascending points at each level
for k = 1:nk
    if na(k) > 0
        bam(k) = (sum(sum(ba(:,:,k),1),2))/na(k);
    end
end

% calculate avg. qh of ascending points at each level
for k = 1:nk
    if na(k) > 0
        qham(k) = (sum(sum(qha(:,:,k),1),2))/na(k);
    end
end



% store cbmf matrix to cbmassfx array
asmf.(fieldname).amf = mf;   % total ascending mass flux at each level
asmf.(fieldname).na = na;   % total number of ascending points at each level
asmf.(fieldname).aat = aat;   % total area of asc. points at each level
asmf.(fieldname).wam = wam;   % avg. vertical velocity of asc. points at each level
asmf.(fieldname).bam = bam;   % avg. buoyancy of asc. points at each level
asmf.(fieldname).qham = qham;   % avg. qh (hydrometeor) of asc. points at each level
asmf.(fieldname).lclxm = lclxm;   % mean LCL over the hour in m
% asmf.(fieldname).zid = zid;
asmf.(fieldname).xh = xh;   % in km
asmf.(fieldname).xhii = xhii;   
asmf.(fieldname).xhfi = xhfi;
asmf.(fieldname).z = z ;   

% save workspace variables 
if saveop == 1
save(fullfile(fname1,filename1),'asmf');
end

end % end seq = 1:12

toc
