
close all
clear all

%%----------------------------------------------------------------------------
% REMARKS:
% cm1 output files are written every 90 s
% This code calculates profile of avg. frac. entrainment rate over
% convective cores using ice-liquid water potential temperature (theta_il)
% The code generates entm array which stores hourly avg. of the 
% quantities mentioned
% formula of theta_il is obtained from eq. (23) of Bryan and Fritsch (2004)
%%-----------------------------------------------------------------------------

tic

% initialize array for storing ent. matrix for each hour
entm_thil = { };

for seq = 6:6
    
% enter range of time (hour) interested:
% ti = 5; tf = 6;
ti = 0+(seq-1)*1; tf = 1+(seq-1)*1;

dt = 1; % time step for reading cm1 output files
dcmt = 90; % time step between each cm1out file written (s)

% directory to read output files (cm1 files) from:
fname = '/home/tang235/projects/ctb-kirshbau/tang235/dx250_mar07';   % beluga
%fname = '/storage/mtang/dx250_mar07';

% % directory for storing ent_m array and name of matrix:
fname1 = '/home/tang235/projects/ctb-kirshbau/tang235/dx250_mar07/matrices';
%fname1 = '/storage/mtang/dx250_mar07/matrices';
filename1 = sprintf('entm_thil_dx250_mar07.mat');

% field name for structure
fieldname = sprintf('f%d%d',ti,tf);

% enter 1 for saving output matrices and 0 for not saving outputs
saveop = 1;

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

cmt = ti*3600/dcmt+1:dt:tf*3600/dcmt+1; %time for cm1out_00000t.nc file 

% INPUT CONSTANTS FROM constants.F
g = 9.81; % m/s^2
p00 = 1.0e5;   % Pa, ref. prs. 
Rd = 287.04;   % J/kg/K, gas constant for dry air
Rv = 461.5;   % J/kg/K, gas constant for water vapor
Cp = 1005.7;   % J/kg/K, spec. heat of dry air at con. prs.
Cpv = 1870;   % J/kg/K, spec. heat of water vapor at con. prs.
Cpl = 4190;   % J/kg/K, spec. heat of liquid water at con. prs.
Cpi = 2106;   % J/kg/K, spec. heat of ice at con. prs.
epsi = Rd/Rv;
lv = 2501000.0;   % J/kg, latent heat of vaporization
ls = 2834000.0;   % J/kg, latent heat of sublimation
t0 = 273.15;   % K, ref. temp.

for s = cmt(1):dt:cmt(end)

%filename = '/aos/home/mtang/cm1out_000241.nc';
filename = sprintf('cm1out_%06d.nc',s);
cmfile = fullfile(fname,filename);

% read data from output files

z = ncread(cmfile,'z');   % z in km
w = ncread(cmfile,'w');      % w in m/s

qv = ncread(cmfile,'qv');    % qv, water vapor in kg/kg
qc = ncread(cmfile,'qc');    % qc, cloud water in kg/kg
qr = ncread(cmfile,'qr');    % qr, rain 
qs = ncread(cmfile,'qs');    % qs, snow
qi = ncread(cmfile,'qi');    % qi in kg/kg
qg = ncread(cmfile,'qg');    % qg, graupel

prs = ncread(cmfile,'prs'); % pressure in Pa
th = ncread(cmfile,'th'); % potential temp. in K

qt = qv + qc + qr + qs + qi + qg; % total mixing ratio

% get size of each direction
ni=size(qc,1);
nj=size(qc,2);
nk=size(qc,3);

% get sum of qc and qi
qcqi = qc + qi;

% Calculate temperature, temp
exn = (prs./p00).^(Rd/Cp);       %exner function
temp = exn.*th;
tempr = temp.*(1+qv./epsi)./(1+qt);     %density temp.

% % Calculate moist static energy
% z3d = zeros(ni,nj,nk);
% for i = 1:ni
%     for j = 1:nj
%         z3d(i,j,:) = z*1000;   % z in m
%     end
% end
% 
% se = Cp*tempr + g*z3d + Lv.*qv;  % moist static energy

% calculate buoyency using theta_r
thr = th.*(1+qv./epsi)./(1+qt);  % theta_r
% thrm = mean(mean(thr,1),2);  % mean theta_r
% thrm3d = repmat(thrm,[ni nj 1]);

b = g*((thr - thr03d)./thr03d);    % buoyency 
% b = g*((thr - thrm3d)./thrm3d);    % buoyency 

% GET w ONTO HALF LEVEL (z)
nkw = size(w,3);
wh = 0.5*(w(:,:,1:nkw-1) + w(:,:,2:nkw));

% set threshold qth
qth = 0.1;  % threshold qc+qi in g/kg  
qth = qth/1000;  % qth in kg/kg

%-----------------------------------------------------
% calculate ice-liquid water potential temperature
% from eq. (23) of Bryan and Fritsch (2004)

% get chi and gamma
chi = (Rd + Rv.*qt)./(Cp + Cpv.*qt);   % eq.(19)
gamma = Rv.*qt./(Cp + Cpv.*qt);   % eq. (20)

% get functional form of Lv and Ls 
Lv = lv + (Cpv-Cpl).*(tempr-t0);   % Lv(T)
Ls = ls + (Cpv-Cpi).*(tempr-t0);   % Ls(T)

% get liquid water and ice water mixing ratio
ql = qc + qr;   % in kg/kg
qii = qs + qi + qg;   % in kg/kg

% ql = qc;   % in kg/kg
% qii = qi;   % in kg/kg

% get each term in eq.(23) and get theta_il
term1 = tempr.*((p00./prs).^chi);
term2 = (1-(ql+qii)./(epsi+qt)).^chi;
term3 = (1-(ql+qii)./qt).^(-gamma);
term4 = exp((-Lv.*ql-Ls.*qii)./((Cp+Cpv.*qt).*tempr));

thil = term1.*term2.*term3.*term4;   % ice-liquid water potential temp. (K)

% done calculating ice-liquid water potential temp.
%-----------------------------------------------------
%-----------------------------------------------------

% calculate frac. entrainment rate using eq.(18) of de Rooy et al. (2013)

% initialize thil_c, thil_e, and fc matrices at the beginning 
if s == cmt(1)
    thilc = zeros(ni,nj,nk);  % th_il of convective core
    thile = zeros(ni,nj,nk);  % th_il of environment
    fc = zeros(ni,nj,nk);  % freq. counter for convective core at each point
end

for i = 1:ni
    for j = 1:nj
        for k = 1:nk
            if qcqi(i,j,k) > qth && wh(i,j,k) > 0 && b(i,j,k) > 0   % conditions to be a convective core
                thilc(i,j,k) = thilc(i,j,k) + thil(i,j,k);
                fc(i,j,k) = fc(i,j,k) + 1; % increment count by 1
            end
        end
    end
end

% calculate thil_e (th_il of environment)
thile = thile + thil;

% sprintf('cm1out_%06d.nc read and hc and he at s = %d calculated \n tih is %d and tfh is %d',s,s,ti,tf)
sprintf(['cm1out_%06d.nc read and at s = %d  \n max. th_il_c = %d, max. th_il_e = %d and tih is %d and tfh is %d'],...
    s,s,max(max(max(thilc))),max(max(max(thile))),ti,tf)

end  % end s = cmt

% calculate avg. thil_c over number convective cores at each level
thilcm = zeros(1,nk);  % profile of thil_c vs z
thilem = zeros(1,nk);  % profile of thil_e vs z
ncc = zeros(1,nk);  % total number of convective cores at each level

% get total number of convective cores at each level 
for k = 1:nk
    ncc(k) = sum(sum(fc(:,:,k),1),2); 
end

% calculcate avg. thil_c over convective core points at each level
for k = 1:nk
    if ncc(k) > 0
        thilcm(k) = (sum(sum(thilc(:,:,k),1),2))/ncc(k);
    end
end

% calculate avg. thil_e at each level
for k = 1:nk
 thilem(k) = (sum(sum(thile(:,:,k),1),2))./(length(cmt)*(ni*nj));
end

% calcualte dthil_c/dz
dthilcmdz = zeros(1,nk);

dz = 1000.*(z(2:end) - z(1:end-1));  % dz in m

% % central difference for thilcm(1,2:nk-1)
% for k = 2:nk-1
%     dthilcmdz(1,k) = (thilcm(1,k+1) - thilcm(1,k-1))./(dz(k) + dz(k-1));
% end
% 
% % forward difference for thilcm(1,1):
% dthilcmdz(1,1) = (thilcm(1,2) - thilcm(1,1))./dz(1);
% 
% % backward difference for thilcm(1,nk):
% dthilcmdz(1,nk) = (thilcm(1,nk) - thilcm(1,nk-1))./dz(end);

%-------------------------------------------------
% get derivatives using fourth-order formulae

% central difference for thilcm(1,2:nk-1)
for k = 3:nk-2
    dthilcmdz(1,k) = (-thilcm(1,k+2) + 8*thilcm(1,k+1) - 8*thilcm(1,k-1) + thilcm(1,k-2))./(12*dz(k));
end

% forward difference for thilcm(1,1:2):
for k = 1:2
    dthilcmdz(1,k) = (-thilcm(1,k+2) + 4*thilcm(1,k+1) - 3*thilcm(1,k))./(2*dz(k));
end

% backward difference for thilcm(1,nk-1:nk):
for k = nk-1:nk
    dthilcmdz(1,k) = (3*thilcm(1,k) - 4*thilcm(1,k-1) + thilcm(1,k-2))./(2*dz(k-1));
end

%-------------------------------------------------

% calculate fractional entrainment rate
ent = -dthilcmdz./(thilcm-thilem);

% assign ent matrix to array named ent 
entm_thil.(fieldname).ent = ent;   % profile of frac. ent. rate of conv. cores 
entm_thil.(fieldname).thilcm = thilcm;   % profile of avg. theta_il over conv. cores 
entm_thil.(fieldname).thilem = thilem;   % profile of avg. theta_il of the environment
entm_thil.(fieldname).z = z;

% %%WRITE entm AS OUTPUT MATRIX FOR FURTHER ANALYSIS
if saveop == 1
    save(fullfile(fname1,filename1),'entm_thil');
end

end  % end for seq = 1:6

toc



