close all
clear all

%%----------------------------------------------------------------------------
% REMARKS:
% If cm1out files are written every 90s, set dcmt=90;
% If cm1out files are written evbery 180s, set dcmt=180;
% The terms are time avg. (hourly) of y-avg. of:
% (ia) u, (ib) w,
% (2) exn-star,
% (3a) subgrid tke, (3b) resolved tke, (3c) total tke in kg/m/s^2 
% terms analyzed would give insight into the relative importance of each
% term in the Dw/Dt eqn. used by dwdt_eqn_1.m
%%-----------------------------------------------------------------------------
tic

% initialize array for storing time avg. of y-avg. of each term in dw/dt eqn for each hour
terms = { };

for seq = 1:12  %((ti,tf) = (0,1)) 1:12 for all hours
    
% enter range of time (hour) interested:
% ti = 0; tf = 12;
ti = 0+(seq-1)*1; tf = 1+(seq-1)*1;

dt = 1; % time step for reading cm1 output files
dcmt = 90; % time step between each cm1out file written (s)

% set range of x- and y-domain intersted:
xmin = -60; % xmin in km
xmax = 60; % xmax in km

ymin = -90; % ymin in km, 5 originnally
ymax = 90; % ymax in km, 5 originally
	
% directory to read output files (cm1 files) from:
fname = '/gs/project/qrw-161-ad/mtang/dx250_mar07';
% fname = '/storage/mtang/dx250_mar07';

% % directory for storing terms array and name of array:
fname1 = '/gs/project/qrw-161-ad/mtang/dx250_mar07/matrices';
% fname1 = '/storage/mtang/dx250_mar07/matrices';
filename1 = sprintf('myuw_exnstar_tke_dx250_mar07.mat');   % terms array

% enter 1 for saving output matrices and 0 for not saving outputs
saveop = 1;

%%---------------------------------------------------------------------------
%%---------------------------------------------------------------------------

% field name for structure 
fieldname = sprintf('f%d%d',ti,tf);

% series for reading cm1 output files
cmt = ti*3600/dcmt+1:dt:tf*3600/dcmt+1; %time for cm1out_00000t.nc file 

% INPUT CONSTANTS FROM constants.F
g = 9.81; % m/s^2
p00 = 1.0e5;
Rd = 287.04;
Cp = 1005.7;
Rv = 461.5;
epsi = Rd/Rv;
Lv = 2501000.0;

%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

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

prs = ncread(cmfile,'prs'); % pressure in Pa
th = ncread(cmfile,'th'); % potential temp. in K

% calculate density using density temperature
exn = (prs./p00).^(Rd/Cp);       %exner function
temp = exn.*th;
tempr = temp.*(1+qv./epsi)./(1+qt);     %density temp.

rho = prs./(Rd*(tempr));       %density

% calculate horizontal mean of initial state rho
r0 = mean(mean(rho,1),2);

% get half-level size
ni=size(qc,1);
nj=size(qc,2);
nk=size(qc,3);

% repmat r_0 into 3-d matrix
r03d = repmat(r0, [ni nj 1]);   % initial state density vs z

% calculate horizontal mean of initial state exn function 
% after reshaping each layer of exn into column vector
exn0 = zeros(1,1,nk);
for k = 1:nk
    exn0(k) = mean(reshape(exn(:,:,k),[ni*nj 1]));  
end

% repmat exn0 into 3-d matrix
exn03d = repmat(exn0, [ni nj 1]);   % initial state hori. avg. (whole domain) exner function vs z

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

%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

% begin reading required files

for s = cmt(1):dt:cmt(end)

filename = sprintf('cm1out_%06d.nc',s);
cmfile = fullfile(fname,filename);

% read data from output files
xf = ncread(cmfile,'xf');   % xf in km
xh = ncread(cmfile,'xh');   % xh in km
yf = ncread(cmfile,'yf');   % yf in km
yh = ncread(cmfile,'yh');   % yh in km
z = ncread(cmfile,'z');   % z in km
zf = ncread(cmfile,'zf');   % zf in km

u = ncread(cmfile,'u');      % u in m/s
v = ncread(cmfile,'v');      % v in m/s
w = ncread(cmfile,'w');      % w in m/s

qv = ncread(cmfile,'qv');    % qv, water vapor in kg/kg
qc = ncread(cmfile,'qc');    % qc, cloud water in kg/kg
qr = ncread(cmfile,'qr');    % qr, rain 
qs = ncread(cmfile,'qs');    % qs, snow
qi = ncread(cmfile,'qi');    % qi, cloud ice in kg/kg
qg = ncread(cmfile,'qg');    % qg, graupel 

qt = qv + qc + qr + qs + qi + qg; % total mixing ratio

prs = ncread(cmfile,'prs'); % pressure in Pa
th = ncread(cmfile,'th'); % potential temp. in K

% get half-level size
ni=size(qc,1);
nj=size(qc,2);
nk=size(qc,3);

% calculate density using density temperature
exn = (prs./p00).^(Rd/Cp);       %exner function
temp = exn.*th;
tempr = temp.*(1+qv./epsi)./(1+qt);     %density temp.

rho = prs./(Rd*(tempr));       %density

% calculate exner function and its perturbation
% exn = (prs./p00).^(Rd/Cp);       %exner function
% exnm = mean(mean(exn,1),2);   % horizontal mean of exn
% exnm3d = repmat(exnm,[ni nj 1]);
% exnstar = exn - exnm3d;   % perturbation exn w.r.t horizontal mean (exn-star)
exnstar = exn - exn03d;   % exn-star: perturbation of exn w.r.t intial state horizontal mean exn 

% calculate buoyency using theta_r
thr = th.*(1+qv./epsi)./(1+qt);  % theta_r
thrm = mean(mean(thr,1),2);  % mean theta_r (horizontal mean)
thrm3d = repmat(thrm,[ni nj 1]);

%b = g*((thr - thrm3d)./thrm3d);    % buoyency 

% get buoyancy using initial state theta_r
b = g*((thr - thr03d)./thr03d);    % buoyency 

% get half-level size
ni=size(qc,1);
nj=size(qc,2);
nk=size(qc,3);

%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

% get u, v, and w onto half levels
uh = 0.5*(u(1:end-1,:,:) + u(2:end,:,:));
vh = 0.5*(v(:,1:end-1,:) + v(:,2:end,:));
wh = 0.5*(w(:,:,1:end-1) + w(:,:,2:end));

% get dxh
dxh = xh(2:end) - xh(1:end-1); % in km

% get indices of xmin and xmax on xh
xhii = max(find(abs(xh-xmin) < dxh(1)/1.5));
xhfi = min(find(abs(xh-xmax) < dxh(1)/1.5));

% reduce x-domain of uh, vh, and wh, into chosen region
uhr = uh(xhii:xhfi,:,:);
vhr = vh(xhii:xhfi,:,:);
whr = wh(xhii:xhfi,:,:);

% get dxh, dyh, and dzh (half-level step size)
dxh = 1000.*(xh(2:end) - xh(1:end-1));   % in m
dyh = 1000.*(yh(2:end) - yh(1:end-1));   % in m
dzh = 1000.*(z(2:end) - z(1:end-1));   % in m, from half-level z

% reduce dxh, dyh, and dzh into chosen region
dxhr = dxh(xhii:xhfi-1);   % in m
dyhr = dyh;   % in m
dzhr = dzh;   % in m

%----------------------------------------------------------------------------

%----------------------------------------------------------------------------

% (VIII) CALCULATE Y-AVG. U, W, and Pi* (exn-star)

% get y-avg. u and w
myuhr = mean(uhr,2); 
mywhr = mean(whr,2); 

% reshape myuhr and mywhr into 2-d (xhfi-xhii+1 * nk) matrices
myuhrr = reshape(myuhr,[xhfi-xhii+1 nk]);   % y-avg. uhr
mywhrr = reshape(mywhr,[xhfi-xhii+1 nk]);   % y-avg. whr

% exn-star: perturb. of exn w.r.t. initial state horizontal mean 
% reduce exnstar into chosen region and
% get y-avg. of exnstar (y-avg. of perturb. of exn w.r.t. horizontal mean)
% i.e. get y-avg. of exn-star
myexnstar = mean(exnstar(xhii:xhfi,:,:),2);

% reshape myexnstar into 2-d (xhfi-xhii+1 * nk) matrix
myexnstarr = reshape(myexnstar,[xhfi-xhii+1 nk]);   % y-avg. uhr

% initialize sum of y-avg u, w, and exn-star matrices at the beginnig of time
if s == cmt(1)
    myuhrs = zeros(xhfi-xhii+1,nk);   % sum of y-avg. u in half-levels and in chosen region
    mywhrs = zeros(xhfi-xhii+1,nk);   % sum of y-avg. w in half-levels and in chosen region
    myexnstars = zeros(xhfi-xhii+1,nk);   % sum of y-avg. exn-star in half-levels and in chosen region
end

% add reshaped y-avg. u to myuhrs matrix and reshaped y-avg. w to mywhrs,
% same for y-avg. exn-star
myuhrs = myuhrs + myuhrr.*dcmt;
mywhrs = mywhrs + mywhrr.*dcmt;
myexnstars = myexnstars + myexnstarr.*dcmt;

% sprintf(['cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d \n,'...
%     'max. of y-avg. u is %d, max. of y-avg. w is %d, max. of y-avg. exn-star is %d,',s,ti,tf,...
%     s,max(max(myuhrs)),max(max(max(mywhrs))),max(max(max(myexnstars)))])

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of y-avg. u is %d',s,ti,tf,s,max(max(myuhrs)))
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of y-avg. w is %d',s,ti,tf,s,max(max(mywhrs)))
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of y-avg. exn-star is %d',s,ti,tf,s,max(max(myexnstars)))

% calculate time avg. of y-avg. u, w, and exn-star at the end of time
if s == cmt(end)
myuhrm = myuhrs./(length(cmt)*dcmt);
mywhrm = mywhrs./(length(cmt)*dcmt);
myexnstarm = myexnstars./(length(cmt)*dcmt);

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of y-avg. u is %d',s,ti,tf,s,max(max(myuhrm)))
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of y-avg. w is %d',s,ti,tf,s,max(max(mywhrm)))
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of y-avg. exn-star is %d',s,ti,tf,s,max(max(myexnstarm)))

end

% END CALCULATING Y-AVG. U, W, AND EXNER-STAR
%----------------------------------------------------------------------------

% (IXA) CALCULATING Y-AVG. TKE (SUBGRID)

% get subgrid tke from output files
tke = ncread(cmfile,'tke'); % subgrid tke in m^2/s^2
sgtke = tke;

% get subgrid tke onto half-level (z)
sgtkeh = 0.5*(sgtke(:,:,1:end-1) + sgtke(:,:,2:end));

% reduce half-level subgrid tke into chosen region
sgtkehr = sgtkeh(xhii:xhfi,:,:);   % in m^2/s^2

% reduce density rho into chosen region
rhor = rho(xhii:xhfi,:,:);

% get subgrid tke in kg/m/s^2 by multiplying subgrid tke with rho
rhosgtker = rhor.*sgtkehr;   

% get y-avg. subgrid tke 
myrhosgtke = mean(rhosgtker,2);   % in kg/m/s^2

% reshape y-avg. subgrid tke in 2-d (xhfi-xhii+1 * nk) matrix
myrhosgtker = reshape(myrhosgtke,[xhfi-xhii+1 nk]);   % in kg/m/s^2

% initialize sum of subgrid tke matrix at the beginning of time
if s == cmt(1)
    myrhosgtkes = zeros(xhfi-xhii+1,nk);
end

% add reshaped y-avg. resolved tke matrix to its sum
myrhosgtkes = myrhosgtkes + myrhosgtker.*dcmt;

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of y-avg. subgrid tke is %d',s,ti,tf,s,max(max(myrhosgtkes)))

% calculate time avg. of y-avg. subgrid tke at the end of time
if s == cmt(end)
myrhosgtkem = myrhosgtkes./(length(cmt)*dcmt);   % in kg/m/s^2

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of y-avg. subgrid tke is %d',s,ti,tf,s,max(max(myrhosgtkem)))

end

% END CALCULATING Y-AVG. SUBGRID TKE
%----------------------------------------------------------------------------

% (IXB) CALCULATING Y-AVG. TKE (RESOLVED)

% y-avg. uhr, vhr, and whr
myuhr = mean(uhr,2);
myvhr = mean(vhr,2); 
mywhr = mean(whr,2);

% repmat y-avg. uhr, vhr, and whr into 3-d matrices
myuhr3d = repmat(myuhr,[1 nj 1]);
myvhr3d = repmat(myvhr,[1 nj 1]);
mywhr3d = repmat(mywhr,[1 nj 1]);

% get perturb. of u from y-avg. u, same for vhr, and whr
pmyuhr = uhr - myuhr3d;   % uh - y-avg. uh
pmyvhr = vhr - myvhr3d;   % vh - y-avg. vh
pmywhr = whr - mywhr3d;   % wh - y-avg. wh

% get resolved tke in chosen region
retke = 0.5.*(pmyuhr.^2 + pmyvhr.^2 + pmywhr.^2);   % in m^2/s^2

% reduce density rho into chosen region
rhor = rho(xhii:xhfi,:,:);

% get resolved TKE in kg/m/s^2 by multiplying resolved tke with rho
rhoretke = rhor.*retke;

% get y-avg. resolved tke (in kg/m/s^2)
myrhoretke = mean(rhoretke,2);

% reshape y-avg. resolved tke in 2-d (xhfi-xhii+1 * nk) matrix
myrhoretker = reshape(myrhoretke,[xhfi-xhii+1 nk]);

% initialize sum of resolved tke matrix at the beginning of time
if s == cmt(1)
    myrhoretkes = zeros(xhfi-xhii+1,nk);
end

% add reshaped y-avg. resolved tke matrix to its sum
myrhoretkes = myrhoretkes + myrhoretker.*dcmt;

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of y-avg. resolved tke is %d',s,ti,tf,s,max(max(myrhoretkes)))

% calculate time avg. of y-avg. resolved tke at the end of time
if s == cmt(end)
myrhoretkem = myrhoretkes./(length(cmt)*dcmt);   % in kg/m/s^2

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of y-avg. resolved tke is %d',s,ti,tf,s,max(max(myrhoretkem)))

end

% END CALCULATING Y-AVG. RESOLVED TKE
%----------------------------------------------------------------------------

% (IXC) CALCULATING Y-AVG. TKE (TOTAL)

% get total tke in kg/m/s^2 (subgrid + resolved) in the chosen region
rhottke = rhosgtker + rhoretke;

% get y-avg. total tke in chosen region
myrhottke = mean(rhottke,2);

% reshape y-avg. total tke in 2-d (xhfi-xhii+1 * nk) matrix
myrhottker = reshape(myrhottke,[xhfi-xhii+1 nk]);

% initialize sum of total tke matrix at the beginning of time
if s == cmt(1)
    myrhottkes = zeros(xhfi-xhii+1,nk);
end

% add reshaped y-avg. total tke matrix to its sum
myrhottkes = myrhottkes + myrhottker.*dcmt;

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of y-avg. total tke is %d',s,ti,tf,s,max(max(myrhottkes)))

% calculate time avg. of y-avg. total tke at the end of time
if s == cmt(end)
myrhottkem = myrhottkes./(length(cmt)*dcmt);   % in kg/m/s^2

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of y-avg. total tke is %d',s,ti,tf,s,max(max(myrhottkem)))

end

% END CALCULATING Y-AVG. TOTAL TKE
%----------------------------------------------------------------------------

% get y-avg. pressure

% reduce density rho into chosen region
prsr = prs(xhii:xhfi,:,:);

% get y-avg.pressure 
myprs = mean(prsr,2);   % in kg/m/s^2

% reshape y-avg. pressure tke in 2-d (xhfi-xhii+1 * nk) matrix
myprsr = reshape(myprs,[xhfi-xhii+1 nk]);   % in Pa

% initialize sum of pressure tke matrix at the beginning of time
if s == cmt(1)
    myprsrs = zeros(xhfi-xhii+1,nk);
end

% add reshaped y-avg. pressure matrix to its sum
myprsrs = myprsrs + myprsr.*dcmt;

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of y-avg. pressure is %d',s,ti,tf,s,max(max(myprsrs)))

% calculate time avg. of y-avg. pressure at the end of time
if s == cmt(end)
myprsrsm = myprsrs./(length(cmt)*dcmt);   % in kg/m/s^2

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of y-avg. pressure is %d',s,ti,tf,s,max(max(myprsrsm)))

end





end  % end s = cmt(1):dt:cmt(end)



% save ouput to terms array
terms.(fieldname).myuhrm = myuhrm;  % time avg. of y-avg. of u
terms.(fieldname).mywhrm = mywhrm;  % time avg. of y-avg. of w
terms.(fieldname).myexnstarm = myexnstarm;  % time avg. of y-avg. of exn-star
terms.(fieldname).myrhosgtkem = myrhosgtkem;  % time avg. of y-avg. of rho*subgrid tke each hour (in kg/m/s^2)
terms.(fieldname).myrhoretkem = myrhoretkem;  % time avg. of y-avg. of rho*resolved tke each hour (in kg/m/s^2)
terms.(fieldname).myrhottkem = myrhottkem;  % time avg. of y-avg. of rho*total tke each hour (in kg/m/s^2)
terms.(fieldname).myprsrsm = myprsrsm;  % time avg. of y-avg. of pressure each hour (in Pa)

terms.(fieldname).xh = xh;  % in km
terms.(fieldname).z = z;  % in km
terms.(fieldname).xhii = xhii;  % in km   % index of xmin on xh
terms.(fieldname).xhfi = xhfi;  % in km   % index of xmax on xh


% save workspace variables 
if saveop == 1
save(fullfile(fname1,filename1),'terms');
end


end  % end seq = 1:12

toc
