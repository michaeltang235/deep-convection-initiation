close all
clear 

%%----------------------------------------------------------------------------
% REMARKS:
% This code computes time avg. of each term in w-budget output of model
% every hour
% rhs = (1/n)*0.5*(rh0+rhs_n) + (rhs1+rhs3+...+rhs_n-1), where n is last
% output file number (0,1,2,...,n-1,n)
% There are possbily some coding bugs that model budget buoyancy to be too
% small than normal
% The diff. btw. wbudget_1.m and wbudget.m are that the sfc. layer of the
% buoaycny field is replaced by the sfc. buoy. calculated by ourselves
% Also, sfc. layer of rdamp field is replaced by extrapolatting the 2nd
% layer down to the surface
%%-----------------------------------------------------------------------------
tic

% initialize array for storing time avg. of y-avg. of each term in dw/dt eqn for each hour
terms = { };

for seq = 6:6  %((ti,tf) = (0,1)) 1:12 for all hours
    
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
fname = '/home/tang235/projects/ctb-kirshbau/tang235/prcl_dx250_dec03';
% fname = '/storage/mtang/prcl_dx250_dec03';

% % directory for storing terms array and name of array:
fname1 = '/home/tang235/projects/ctb-kirshbau/tang235/prcl_dx250_dec03/matrices';
% fname1 = '/storage/mtang/prcl_dx250_dec03/matrices';
filename1 = sprintf('wbudget_1_prcl_dx250_dec03.mat');   % terms array

% enter 1 for saving output matrices and 0 for not saving outputs
saveop = 1;

%%---------------------------------------------------------------------------
%%---------------------------------------------------------------------------

%----------------------------------------------------------------------------
% Calculate initial state density potential temperature

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

% series for reading cm1 output files
cmt = ti*3600/dcmt+1:dt:tf*3600/dcmt+1; %time for cm1out_00000t.nc file 

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

qv = ncread(cmfile,'qv');    % qv, water vapor in kg/kg
qc = ncread(cmfile,'qc');    % qc, cloud water in kg/kg
qr = ncread(cmfile,'qr');    % qr, rain 
qs = ncread(cmfile,'qs');    % qs, snow
qi = ncread(cmfile,'qi');    % qi, cloud ice in kg/kg
qg = ncread(cmfile,'qg');    % qg, graupel 

qt = qv + qc + qr + qs + qi + qg; % total mixing ratio

w = ncread(cmfile,'w');      % w in m/s

th = ncread(cmfile,'th'); % potential temp. in K

wb_hadv = ncread(cmfile,'wb_hadv');   % wb_hadv in m/s/s
wb_vadv = ncread(cmfile,'wb_vadv');   % wb_vadv in m/s/s
wb_hidiff = ncread(cmfile,'wb_hidiff');   % wb_hidiff in m/s/s
wb_vidiff = ncread(cmfile,'wb_vidiff');   % wb_hidiff in m/s/s
wb_hturb = ncread(cmfile,'wb_hturb');   % wb_hturb in m/s/s
wb_vturb = ncread(cmfile,'wb_vturb');   % wb_vturb in m/s/s
wb_pgrad = ncread(cmfile,'wb_pgrad');   % wb_pgrad in m/s/s
wb_rdamp = ncread(cmfile,'wb_rdamp');   % wb_rdamp in m/s/s
wb_buoy = ncread(cmfile,'wb_buoy');   % wb_rdamp in m/s/s

% get full-level sizes
ni = size(wb_hadv,1);
nj = size(wb_hadv,2);
nkf = size(wb_hadv,3);

%-----------------------------------------------------
%----------------------------------------------------------------------------
% calculate buoyancy field ourselves and extrapolate it down to surface
% (z=0)

% calculate density potential temperature at current time step
thr = th.*(1+qv./epsi)./(1+qt);  % theta_r

% get buoyancy using density potential temperature
b = g*((thr - thr03d)./thr03d);    % buoyancy

%------------------------------
% get horizontal avg. of buoyancy in the whole domain
% after reshaping each layer of buoy. into column vector
bhm = zeros(1,1,nk);
for k = 1:nk
    bhm(k) = mean(reshape(b(:,:,k),[ni*nj 1]));
end

% repmat hori. avg. (whole domain) of buoy. into [ni* 1 * nk] matrix
bhm3d = repmat(bhm,[ni 1 1]);

% subtract reduced hori. avg. of buoy. from y-avg. buoy. in the chosen region
bdhm = b - bhm3d;   % to remove sheet-like appearance 

%------------------------------

% since dimension of wbudget ouput is (ni,nj,nkf),
% convert b to full level in z, and extrapolate it down to the sfc. (z=0)
% so, get full-level (in zf) buoyancy 
% bf = 0.5.*(bdhm(:,:,1:nk-1) + bdhm(:,:,2:nk));
bf = 0.5.*(b(:,:,1:nk-1) + b(:,:,2:nk));

% get slope between the 1st and 2nd layers of bf (full-level buoyancy)
slope = 0.1./(bf(:,:,2) - bf(:,:,1));   % change in b over 0.1 km

% extrapolate buoyancy field (full-level) down to the surface
bsfc = bf(:,:,1) - 0.1./slope;

% end calculating buoyancy field ourselves
%----------------------------------------------------------------------------

% (IA & IB) CALCULATE y-avg. hori. and vert. adv. term

% get y-avg. hadv and vadv
myhadv = mean(wb_hadv,2);   % y-avg. hori. adv. term
myvadv = mean(wb_vadv,2);   % y-avg. vert. adv. term

% squeeze myhadv and myvadv into ni*nk matrices
myhadv2d = squeeze(myhadv);  
myvadv2d = squeeze(myvadv);  

% initialize hori. adv. matrix (hori_adv and ver_adv term) at the beginnig of time
if s == cmt(1)
    hadv = zeros(ni,nkf);   % hori. adv. term
    vadv = zeros(ni,nkf);   % vert. adv. term
end

% add 2d y-avg. hori. adv. to hadv, do the same for vadv
hadv = hadv + (myhadv2d).*dcmt;
vadv = vadv + (myvadv2d).*dcmt;

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. hori. adv. is %d',s,ti,tf,s,max(max(hadv)))
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. vert. adv. is %d',s,ti,tf,s,max(max(vadv)))

% calculate time avg. of the adv term at the end of time
if s == cmt(end)
hadvm = hadv./(length(cmt)*dcmt);   % time avg. of y-avg. hori. adv. term
vadvm = vadv./(length(cmt)*dcmt);   % time avg. of y-avg. vert. adv. term

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of hori. adv. term is %d',s,ti,tf,s,max(max(hadvm)))
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of vert. adv. term is %d',s,ti,tf,s,max(max(vadvm)))

end

% rhs = wb_hadv+wb_vadv+wb_hidiff+wb_vidiff+wb_hturb+wb_vturb+wb_pgrad+wb_rdamp+wb_buoy;

% END CALCULATING ADVECTION TERM (hori. and vert.)
%----------------------------------------------------------------------------

% (IIA & IIB) CALCULATE y-avg. hori. and vert. diff. term

% get y-avg. hidiff and vidiff
myhidiff = mean(wb_hidiff,2);   % y-avg. hori. diff. term
myvidiff = mean(wb_vidiff,2);   % y-avg. vert. diff. term

% squeeze myhidiff and myvidiff into ni*nk matrices
myhidiff2d = squeeze(myhidiff);
myvidiff2d = squeeze(myvidiff);  

% initialize hori. diff. matrixat the beginnig of time
if s == cmt(1)
    hidiff = zeros(ni,nkf);   % hori. diff. term
    vidiff = zeros(ni,nkf);   % vert. diff. term
end

% add 2d y-avg. hori. diff. to hidiff, do the same for vidiff
hidiff = hidiff + (myhidiff2d).*dcmt;
vidiff = vidiff + (myvidiff2d).*dcmt;

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. hori. diff. is %d',s,ti,tf,s,max(max(hidiff)))
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. vert. diff. is %d',s,ti,tf,s,max(max(vidiff)))

% calculate time avg. of the diff. term at the end of time
if s == cmt(end)
hidiffm = hidiff./(length(cmt)*dcmt);   % time avg. of y-avg. hori. diff. term
vidiffm = vidiff./(length(cmt)*dcmt);   % time avg. of y-avg. vert. diff. term

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of hori. diff. term is %d',s,ti,tf,s,max(max(hidiffm)))
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of vert. diff. term is %d',s,ti,tf,s,max(max(vidiffm)))

end

% END CALCULATING DIFFUSION TERM (hori. and vert.)
%----------------------------------------------------------------------------

% (IIIA & IIIB) CALCULATE y-avg. hori. and vert. turb. term

% get y-avg. hturb and vturb
myhturb = mean(wb_hturb,2);   % y-avg. hori. turb. term
myvturb = mean(wb_vturb,2);   % y-avg. vert. turb. term

% squeeze myhturb and myvturb into ni*nk matrices
myhturb2d = squeeze(myhturb);
myvturb2d = squeeze(myvturb);  

% initialize hori. turb. matrixat the beginnig of time
if s == cmt(1)
    hturb = zeros(ni,nkf);   % hori. turb. term
    vturb = zeros(ni,nkf);   % vert. turb. term
end

% add 2d y-avg. hori. turb. to hturb, do the same for vturb
hturb = hturb + (myhturb2d).*dcmt;
vturb = vturb + (myvturb2d).*dcmt;

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. hori. turb. is %d',s,ti,tf,s,max(max(hturb)))
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. vert. turb. is %d',s,ti,tf,s,max(max(vturb)))

% calculate time avg. of the hori. and vert. turb. term at the end of time
if s == cmt(end)
hturbm = hturb./(length(cmt)*dcmt);   % time avg. of y-avg. hori. turb. term
vturbm = vturb./(length(cmt)*dcmt);   % time avg. of y-avg. vert. turb. term

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of hori. turb. term is %d',s,ti,tf,s,max(max(hturbm)))
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of vert. turb. term is %d',s,ti,tf,s,max(max(vturbm)))

end

% END CALCULATING TURBULENCE TERM (hori. and vert.)
%----------------------------------------------------------------------------

% (IV) CALCULATE y-avg. pgf term

% get y-avg. pgf
mypgrad = mean(wb_pgrad,2);   % y-avg. pgf term

% squeeze mypgrad into ni*nk matrices
mypgrad2d = squeeze(mypgrad);

% initialize pgrad matrixat the beginnig of time
if s == cmt(1)
    pgrad = zeros(ni,nkf);   % pgf term
end

% add 2d y-avg. pgf to pgrad
pgrad = pgrad + (mypgrad2d).*dcmt;

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. pgf term is %d',s,ti,tf,s,max(max(pgrad)))

% calculate time avg. of the pgf term at the end of time
if s == cmt(end)
pgradm = pgrad./(length(cmt)*dcmt);   % time avg. of y-avg. pgf term

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of pgf term is %d',s,ti,tf,s,max(max(pgradm)))

end

% END CALCULATING PGF TERM
%----------------------------------------------------------------------------

% (V) CALCULATE y-avg. rdamp term

%--------------------------------------------
% extrapolate damping term to the surface
rdslope = 0.1./(wb_rdamp(:,:,3) - wb_rdamp(:,:,2));
wb_rdamp_1 = wb_rdamp;
wb_rdamp_1(:,:,1) = wb_rdamp(:,:,2) - 0.1./rdslope;
%--------------------------------------------

% get y-avg. rdamp
myrdamp = mean(wb_rdamp_1,2);   % y-avg. rdamp term
% myrdamp = mean(wb_rdamp,2);   % y-avg. rdamp term

% squeeze myrdamp into ni*nk matrices
myrdamp2d = squeeze(myrdamp);

% initialize rdamp matrixat the beginnig of time
if s == cmt(1)
    rdamp = zeros(ni,nkf);   % rdamp term
end

% add 2d y-avg. rdamp to rdamp
rdamp = rdamp + (myrdamp2d).*dcmt;

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. rdamp term is %d',s,ti,tf,s,max(max(rdamp)))

% calculate time avg. of the rdamp term at the end of time
if s == cmt(end)
rdampm = rdamp./(length(cmt)*dcmt);   % time avg. of y-avg. rdamp term

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of rdamp term is %d',s,ti,tf,s,max(max(rdampm)))

end

% END CALCULATING RDAMP TERM
%----------------------------------------------------------------------------

% (VI) CALCULATE y-avg. buoy. term

% %---------------------------------------------------
% % add the extrapolated b_sfc to buoyancy field 
wb_buoy_1 = wb_buoy;
wb_buoy_1(:,:,1) = bsfc;
% %---------------------------------------------------

% get y-avg. buoy.
mybuoy = mean(wb_buoy_1,2);   % y-avg. buoy. term

% squeeze mybuoy into ni*nk matrices
mybuoy2d = squeeze(mybuoy);

% initialize buoy matrixat the beginnig of time
if s == cmt(1)
    buoy = zeros(ni,nkf);   % buoy. term
end

% add 2d y-avg. buoy to buoy
buoy = buoy + (mybuoy2d).*dcmt;

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. buoy. term is %d',s,ti,tf,s,max(max(buoy)))

% calculate time avg. of the buoy term at the end of time
if s == cmt(end)
buoym = buoy./(length(cmt)*dcmt);   % time avg. of y-avg. buoy term

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of buoy. term is %d',s,ti,tf,s,max(max(buoym)))

end

% END CALCULATING BUOY. TERM
% %----------------------------------------------------------------------------

% (VII) CALCULATE y-avg. r.h.s.

% get sum of r.h.s at current time step
rhs = wb_hadv + wb_vadv + wb_hidiff + wb_vidiff + wb_hturb + wb_vturb + wb_pgrad + wb_rdamp_1 + wb_buoy_1;
% rhs = wb_hadv + wb_vadv + wb_hidiff + wb_vidiff + wb_hturb + wb_vturb + wb_pgrad + wb_rdamp + wb_buoy;

% get y-avg. r.h.s. at current time step
myrhs = mean(rhs,2);

% squeeze myrhs into ni*nk matrices
myrhs2d = squeeze(myrhs);

% get y-avg. r.h.s. at the first time step
if s == cmt(1)
    myrhsi = myrhs2d;
end

% initialize sum of r.h.s. matrix
if s == cmt(2)
    myrhscrs = zeros(ni,nkf);
end

% get sum of y-avg. r.h.s. of all time steps (except the first and last
% one)
if s ~= cmt(1) && s ~= cmt(end)
    myrhscrs = myrhscrs + myrhs2d;
    
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. y-avg. r.hs. terms is %d',s,ti,tf,s,max(max(myrhscrs)))
    
end

% get y-avg. r.h.s. at the last time step
if s == cmt(end)
    myrhsf = myrhs2d;
end

% get time avg. of y-avg. r.h.s. at the end of time
if s == cmt(end)
myrhsm = (1/(length(cmt)-1)).*(0.5.*(myrhsi + myrhsf) + myrhscrs);

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of r.h.s. terms is %d',s,ti,tf,s,max(max(myrhsm)))
end

% END CALCULATING R.H.S. OF EQN.
%----------------------------------------------------------------------------
% % % ----------------------------------------------------------------------------

% (VIII) CALCULATE d/dt(y-avg. w)

% get y-avg. w at initial time (from reduced w in chosen region)
if s == cmt(1)
mywi = mean(w,2);

% squeeze y-avg. w into ni*nkf matrix
mywir = squeeze(mywi);

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of y-avg. w at initial time is %d',s,ti,tf,s,max(max(mywir)))

end

% get y-avg. w at final time (from reduced w in chosen region)
if s == cmt(end)
mywf = mean(w,2);

% squeeze y-avg. w into ni*nk matrix
mywfr = squeeze(mywf);

% get d/dt(y-avg. w) at the end of time
ddtwm = (mywfr - mywir)./((cmt(end)-cmt(1))*dcmt);

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of y-avg. w at final time is %d',s,ti,tf,s,max(max(mywfr)))

end
% 
% % END CALCULATING d/dt(Y-AVG. w)
% % %----------------------------------------------------------------------------
% % %----------------------------------------------------------------------------
% 
% % (VIII) CALCULATE d/dt(Y-AVG. w) VERSION 2
% 
% % get y-avg. w
% myw = mean(w,2);
% 
% % squeeze mywhr into 2-d (ni*nkf) matrix
% mywr = squeeze(myw);
% 
% sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of y-avg. w is %d',s,ti,tf,s,max(max(mywr)))
% 
% % initialize mywhrar array at the beginning of time to store y-avg. wh at all time steps
% if s == cmt(1)
%     mywhrar = zeros(ni,nkf,length(cmt));
% end
% 
% % get position of current time step in array (k^th dimension)
% ark = find(s==cmt);
% 
% % update mywhrar array with y-avg. wh at current time step
% mywhrar(:,:,ark) = mywr;
% 
% % get dw/dt array at each time step at the end of time
% if s == cmt(end)
%     ddtwar = zeros(ni,nkf,length(cmt));
%     
% % get dw/dt at each time step using second order formulae
% ddtwar(:,:,1) = (mywhrar(:,:,2) - mywhrar(:,:,1))./(dcmt*dt);
% for k = 2:length(cmt)-1
%     ddtwar(:,:,k) = (mywhrar(:,:,k+1) - mywhrar(:,:,k-1))./(2*dcmt*dt);
% end
% ddtwar(:,:,end) = (mywhrar(:,:,end) - mywhrar(:,:,end-1))./(dcmt*dt);
% 
% % get time avg. of dwdt at the end of time
% ddtwm = (sum(ddtwar,3).*dcmt)./(length(cmt)*dcmt);
% 
% sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of dwdt is %d',s,ti,tf,s,max(max(ddtwm)))
% 
% end   % end if s == cmt(end)
% 
% % END CALCULATING d/dt(Y-AVG. w)
% % %----------------------------------------------------------------------------


end   % end for s = cmt

%------------------------------
% save ouput to terms array
terms.(fieldname).hadvm = hadvm;  % time avg. of y-avg. of hori. adv. each hour
terms.(fieldname).vadvm = vadvm;  % time avg. of y-avg. of vert. adv. each hour

terms.(fieldname).hidiffm = hidiffm;  % time avg. of y-avg. of hori. diff. each hour
terms.(fieldname).vidiffm = vidiffm;  % time avg. of y-avg. of vert. diff. each hour

terms.(fieldname).hturbm = hturbm;  % time avg. of y-avg. of hori. turb. each hour
terms.(fieldname).vturbm = vturbm;  % time avg. of y-avg. of vert. turb. each hour

terms.(fieldname).pgradm = pgradm;  % time avg. of y-avg. of pgrad each hour
terms.(fieldname).rdampm = rdampm;  % time avg. of y-avg. of rdamp each hour
terms.(fieldname).buoym = buoym;  % time avg. of y-avg. of buoy. each hour

terms.(fieldname).myrhsm = myrhsm;  % time avg. of y-avg. of r.h.s. terms each hour
terms.(fieldname).ddtwm = ddtwm;  % time avg. of y-avg. of dw/dt each hour

terms.(fieldname).xh = xh;  % in km
terms.(fieldname).zf = zf;  % in km


% save workspace variables 
if saveop == 1
save(fullfile(fname1,filename1),'terms');
end


end   % end seq = 6:6
