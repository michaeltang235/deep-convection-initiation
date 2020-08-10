close all
clear 

%%----------------------------------------------------------------------------
% REMARKS:
% This code computes hourly time avg. of:
% (ia) hori. turb. (Th, subgrid) from w-budget 
% (ib) vert. turb. (Tv, subgrid) from w-budget
% (ii) Dw (subgrid) turb. flux from cm1 eqn. doc. (n.a. in this version)
% (iii) resolved turb. flux of w from own calculation
% get total turb. flux of w in the following ways
% (iv) Th+Tv-resolved flux of w
% (v) Dw-resolved flux of w (n.a. in this version)
% (vi) y-avg. w

% The diff. btw. turbflux_w_2.m and turbflux_w.m is that it considers the
% initial state density (rho_0) when calculating the resolved
% turb. flux (only the vertical) in the w-momentum eq.
% i.e. formula for resolved turb. flux of w' is different from
% turbflux_w_1.m
%%-----------------------------------------------------------------------------
tic

% initialize array for storing time avg. of y-avg. of each term in dw/dt eqn for each hour
terms = { };


for seq = 6:6  %((ti,tf) = (0,1)) 1:12 for all hours
    
% enter range of time (hour) interested:
% ti = 0; tf = 12;
ti = 0+(seq-1)*1; tf = 1+(seq-1)*1;

dt = 4; % time step for reading cm1 output files
dcmt = 90; % time step between each cm1out file written (s)

% set range of x- and y-domain intersted:
xmin = -10; % xmin in km
xmax = 10; % xmax in km
	
% directory to read output files (cm1 files) from:
fname = '/home/tang235/projects/ctb-kirshbau/tang235/prcl_dx250_dec03';   % beluga
% fname = '/storage/mtang/prcl_dx250_dec03';

% % directory for storing terms array and name of array:
fname1 = '/home/tang235/projects/ctb-kirshbau/tang235/prcl_dx250_dec03/matrices';
% fname1 = '/storage/mtang/prcl_dx250_dec03/matrices';
filename1 = sprintf('turbflux_w_2_prcl_dx250_dec03.mat');   % terms array

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

% PRELIMINARY: GET INITIAL STATE DENSITY

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
ni=size(th,1);
nj=size(th,2);
nk=size(th,3);

% repmat r_0 into 3-d matrix
r03d = repmat(r0, [ni nj 1]);   % initial state density vs z

% DONE GEETING 3-D INITIAL STATE DENSITY 
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

% wb_hadv = ncread(cmfile,'wb_hadv');   % wb_hadv in m/s/s
% wb_vadv = ncread(cmfile,'wb_vadv');   % wb_vadv in m/s/s
% wb_hidiff = ncread(cmfile,'wb_hidiff');   % wb_hidiff in m/s/s
% wb_vidiff = ncread(cmfile,'wb_vidiff');   % wb_hidiff in m/s/s
wb_hturb = ncread(cmfile,'wb_hturb');   % wb_hturb in m/s/s
wb_vturb = ncread(cmfile,'wb_vturb');   % wb_vturb in m/s/s
% wb_pgrad = ncread(cmfile,'wb_pgrad');   % wb_pgrad in m/s/s
% wb_rdamp = ncread(cmfile,'wb_rdamp');   % wb_rdamp in m/s/s
% wb_buoy = ncread(cmfile,'wb_buoy');   % wb_rdamp in m/s/s

% get full-level size on the vertical
nkf = size(w,3);

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
% %--------------------------------------------------------------------------
% (IA & IB) GET SUBGRID TURB. FLUX (Th and Tv) FROM W-BUDGET OUTPUT: 

% CALCULATE y-avg. hori. and vert. turb. term

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


% get hturbm and vturbm into chosen region
hturbmr = hturbm(xhii:xhfi,:);   % hturbm in chosen region
vturbmr = vturbm(xhii:xhfi,:);   % vturbm in chosen region

% get sum of subgrid turb. flux (Th+Tv) of w
thtvmr = hturbmr + vturbmr;

% get subgrid turb. flux of w (Th, Tv, and Th+Tv) from w-budget into half-levels
hturbmrr = 0.5.*(hturbmr(:,1:end-1) + hturbmr(:,2:end));   % get Th into half-level
vturbmrr = 0.5.*(vturbmr(:,1:end-1) + vturbmr(:,2:end));   % get Tv into half-level
thtvmrr = 0.5.*(thtvmr(:,1:end-1) + thtvmr(:,2:end));   % get Th+Tv into half-level

end   % end if s == cmt(end)

% END GETTING SUBGIRD TURB. FLUX OF W (Th+Tv)
% %--------------------------------------------------------------------------
% 
% % (II) CALCULATE Y-AVERGAE D_W: subgrid flux (own calculation)
% 
% % reduce initial state density (r_0) into chosen region
% r0r = r03d(xhii:xhfi,:,:);
% 
% % get du/dz, dv/dz, dw/dx, dw/dy, and dw/dz using derivative_1.m
% % (second-order formulae)
% dudz = derivative_1(uhr,dzhr,3,1);
% 
% dvdz = derivative_1(vhr,dzhr,3,1);
% 
% dwdx = derivative_1(whr,dxhr,1,1);
% dwdy = derivative_1(whr,dyhr,2,1);
% dwdz = derivative_1(whr,dzhr,3,1);
% 
% % get vertical eddy mixing coefficient, Kmv, and Kmh
% kmv = ncread(cmfile,'kmv');   % kmv in m^2/s
% kmh = ncread(cmfile,'kmh');   % kmh in m^2/s
% 
% % bring kmv onto half-level and reduce kmv into chosen region
% kmvh = 0.5.*(kmv(:,:,1:end-1) + kmv(:,:,2:end));
% kmvhr = kmvh(xhii:xhfi,:,:);
% 
% % get tau13, tau23, and tau33
% tau13 = kmvhr.*(dudz + dwdx);
% tau23 = kmvhr.*(dvdz + dwdy);
% tau33 = kmvhr.*(dwdz + dwdz);
% 
% % get r_0*tau,
% rt13 = r0r.*tau13;
% rt23 = r0r.*tau23;
% rt33 = r0r.*tau33;
% 
% % get d/dx(r_0*tau13), d/dy(r_0*tau23), and d/dz(r_0*tau33) using
% % derivative_1.m (second-order formulae)
% ddxrtau13 = derivative_1(rt13,dxhr,1,1);
% ddyrtau23 = derivative_1(rt23,dyhr,2,1);
% ddzrtau33 = derivative_1(rt33,dzhr,3,1);
% 
% % get Dw
% dw = (1/r0r).*(ddxrtau13 + ddyrtau23 + ddzrtau33);
% 
% % get y-avg. of Dw and reshape it into xhfi-xhii+1*nk matrix
% mydw = mean(dw,2);   % y-avg. of dw
% mydwr = reshape(mydw,[xhfi-xhii+1 nk]);
% 
% % initialize subgrid flux matrix at the beginning of time
% if s == cmt(1)
%     sgfx = zeros(xhfi-xhii+1,nk);   % subgrid flux due to Dw term
% end
% 
% % add reshaped y-avg. of Dw to subgrid flux matrix
% sgfx = sgfx + mydwr.*dcmt;
% 
% sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of subgrid flux is %d',s,ti,tf,s,max(max(sgfx)))
% 
% % calculate time avg. of the subgrid flux term at the end of time
% if s == cmt(end)
% sgfxm = sgfx./(length(cmt)*dcmt);
% 
% sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of subgrid flux is %d',s,ti,tf,s,max(max(sgfxm)))
% end
% 
% % END CALCULATING Y-AVERGAE D_W
% %----------------------------------------------------------------------------

% (III) RESOLVED TURB. FLUX OF W:

% CALCULATE DIV . y-avg.(VEl'w'): turb. flux

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

% get perturb. uh * perturb. wh, perturb. vh * perturb. wh, and 
% perturb. wh * perturb. wh
pmyuhrpmywhr = pmyuhr.*pmywhr;   % perturb. of uh * perturb. wh
pmyvhrpmywhr = pmyvhr.*pmywhr;   % perturb. of vh * perturb. wh
pmywhrpmywhr = pmywhr.*pmywhr;   % perturb. of wh * perturb. wh

% get rho_0 into chosen region
r03dr = r03d(xhii:xhfi,:,:);   % reduced initial state rho (rho_0) in chosen region

%--------------------
% since rho_0 varies in z, add it to the z-comp. of res. turb. flux of w'
rho0pmywhrpmywhr = r03dr.*pmywhrpmywhr;
%--------------------

% get d/dx(y-avg. of perturb. of uh * perturb. wh)
% get d/dy(y-avg. of perturb. of vh * perturb. wh)
% get d/dz[y-avg. of (rho0 * perturb. of wh * perturb. wh)], 
% using derivative_1.m (second-order formulae)

ddxpmyuhrpmywhr = derivative_1(pmyuhrpmywhr,dxhr,1,1);
ddypmyvhrpmywhr = derivative_1(pmyvhrpmywhr,dyhr,2,1);
ddzrho0pmywhrpmywhr = derivative_1(rho0pmywhrpmywhr,dzhr,3,1);   % rho_0 exits in z-comp. 

% get y-avg. of d/dx(perturb. of uh * perturb. wh), same for vh and wh
myddxpmyuhrpmywhr = mean(ddxpmyuhrpmywhr,2);   % y-avg. of d/dx(perturb. of uh * perturb. wh)
myddypmyvhrpmywhr = mean(ddypmyvhrpmywhr,2);   % y-avg. of d/dy(perturb. of vh * perturb. wh)
myrho0ddzrho0pmywhrpmywhr = mean(((1/r03dr).*ddzrho0pmywhrpmywhr),2);   % y-avg. of [(1/rho_0) * d/dz(rho_0 * perturb. of wh * perturb. wh)]

% squeeze the (d/dx, d/dy, and 1/rho0*d/dz(...) matrices to xhfi-xhii+1*nk matrices
myddxpmyuhrpmywhrr = squeeze(myddxpmyuhrpmywhr);
myddypmyvhrpmywhrr = squeeze(myddypmyvhrpmywhr);
myrho0ddzrho0pmywhrpmywhrr = squeeze(myrho0ddzrho0pmywhrpmywhr);  

% initialize resolevd turb. flux matrix at the beginning of time
if s == cmt(1)
    turbfxresh = zeros(xhfi-xhii+1,nk);   % hori. resolved turb. flux, contains factor of (-1)
    turbfxresv = zeros(xhfi-xhii+1,nk);   % vert. resolved turb. flux 
    turbfxres = zeros(xhfi-xhii+1,nk);   % total (hori. + vert.) resolved turb. flux 
end

% add reshaped d/dx, d/dy, and d/dz(...) matrices to turb. flux matrix
turbfxresh = turbfxresh + (-1).*(myddxpmyuhrpmywhrr + myddypmyvhrpmywhrr).*dcmt;
turbfxresv = turbfxresv + (-1).*myrho0ddzrho0pmywhrpmywhrr.*dcmt;
turbfxres = turbfxres + (-1).*(myddxpmyuhrpmywhrr + myddypmyvhrpmywhrr + myrho0ddzrho0pmywhrpmywhrr).*dcmt;

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of hori. comp. of res. turb. flux of w is %d',s,ti,tf,s,max(max(turbfxresh)))
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of vert. comp. of res. turb. flux of w is %d',s,ti,tf,s,max(max(turbfxresv)))
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of res. turb. flux (hori.vert.) of w is %d',s,ti,tf,s,max(max(turbfxres)))

% calculate time avg. of the resolved turb. flux of w' at the end of time
if s == cmt(end)
turbfxreshm = turbfxresh./(length(cmt)*dcmt);   % time avg. of y-avg. of hori. res. turb. flux of w'
turbfxresvm = turbfxresv./(length(cmt)*dcmt);   % time avg. of y-avg. of vert. res. turb. flux of w'
turbfxresm = turbfxres./(length(cmt)*dcmt);   % time avg. of y-avg. of res. turb. flux of w'

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of hori. res. turb. flux is %d',s,ti,tf,s,max(max(turbfxreshm)))
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of vert. res. turb. flux is %d',s,ti,tf,s,max(max(turbfxresvm)))
sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of res. turb. flux is %d',s,ti,tf,s,max(max(turbfxresm)))
end

% END CALCULATING Y-AVG. of [(-1/rho0) * DIV . (rho0*VEl'w')]: res. turb.
% flux of w'
%----------------------------------------------------------------------------

% (V) CALCULATE TOTAL TURB. FLUX OF W USING (Th+Tv+Resolved TURB. FLUX)

% get time avg. of y-avg. total turb. flux of w by summing 
% time avg. of y-avg. of each term (Th+Tv+resolved turb. flux of w)
if s == cmt(end)
thtvturbfxresm = thtvmrr + turbfxresm;   % Tv+Tv+resolved turb. flux of w  (note: resolved term already has (-1) in it)
end

% END GETTING TOTAL TURB. FLUX OF W USING Th+Tv+resolved turb. flux of w
% %----------------------------------------------------------------------------

% % (VI) CALCULATE TOTAL TURB. FLUX OF W USING (Dw-Resolved TURB. FLUX)
% 
% % get time avg. of y-avg. total turb. flux of w by summing 
% % time avg. of y-avg. of each term (Dw-resolved turb. flux of w)
% if s == cmt(end)
% dwturbfxm = sgfxm - turbfxm;   % Dw-resolved turb. flux of w
% end
% 
% % END GETTING TOTAL TURB. FLUX OF W USING Dw-resolved turb. flux of w
% %----------------------------------------------------------------------------

% (VI) CALCULATE Y-AVG. W

% squeeze mywh into 2-d matrices
mywhr2d = squeeze(mywhr);

% initialize sum of y-avg w matrices at the beginnig of time
if s == cmt(1)
    mywhrs = zeros(xhfi-xhii+1,nk);   % sum of y-avg. w in half-levels and in chosen region
end

% add reshaped y-avg. w to mywhrs
mywhrs = mywhrs + mywhr2d.*dcmt;

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of y-avg. w is %d',s,ti,tf,s,max(max(mywhrs)))

% calculate time avg. of w at the end of time
if s == cmt(end)
mywhrm = mywhrs./(length(cmt)*dcmt);

sprintf('cm1out_%06d.nc read,ti = %d, tf = %d \n at cmt = %d , max. of time avg. of y-avg. w is %d',s,ti,tf,s,max(max(mywhrm)))

end

% END GETTING Y-AVG. W
%-----------------------------------------------------------------------------

end   % end for s = cmt



% %------------------------------

% % save ouput to terms array (matrices are reduced into chosen region)

terms.(fieldname).hturbmrr = hturbmrr;  % time avg. of y-avg. of hori. turb. (subgrid flux from w-budget) each hour
terms.(fieldname).vturbmrr = vturbmrr;  % time avg. of y-avg. of vert. turb. (subgrid flux from w-budget) each hour
terms.(fieldname).thtvmrr = thtvmrr;  % time avg. of y-avg. of Th+Tv (subgrid flux from w-budget) each hour
% terms.(fieldname).sgfxm = sgfxm;  % time avg. of y-avg. of Dw (subgrid flux from cm1 eqn.doc.) each hour

terms.(fieldname).turbfxresm = turbfxresm;  % time avg. of y-avg. of resolved (hori.+vert.) turb. flux of w each hour
terms.(fieldname).turbfxreshm = turbfxreshm;   % time avg. of y-avg. of hori. comp. of resolved turb. flux of w) each hour
terms.(fieldname).turbfxresvm = turbfxresvm;   % time avg. of y-avg. of vert. comp. of resolved turb. flux of w) each hour

% OUTPUT TOTAL TURB. FLUX OF W WHICH WAS OBTAINED USING TWO DIFFERENT WAYS
terms.(fieldname).thtvturbfxresm = thtvturbfxresm;  % time avg. of y-avg. of (Th+Tv+resolved flux of w) each hour
% terms.(fieldname).dwturbfxm = dwturbfxm;  % time avg. of y-avg. of (Dw-resolved flux of w) each hour

terms.(fieldname).mywhrm = mywhrm;  % time avg. of y-avg. of w

terms.(fieldname).xh = xh;  % in km
terms.(fieldname).xhii = xhii;  % index of xmin on xh
terms.(fieldname).xhfi = xhfi;  % index of xmax on xh
terms.(fieldname).z = z;  % in km
terms.(fieldname).zf = zf;  % in km

% save workspace variables 
if saveop == 1
save(fullfile(fname1,filename1),'terms');
end


end   % end seq = 6:6

toc