close all
clear all

%---------------------------------------------------------------
% This code inputs frac. mass flux (fmsfx) and frac. entrainment 
% (entm) matrices and calculates frac. detrainment rate (detm)
% frac. ent. rate is obtained using theta_il.m which calculates the
% quantity using eq.(23) of Bryan and Fritsch (2004)
%---------------------------------------------------------------

% input files

% directory for fmsfx and entm .mat files
% fname = 'aos/home/mtang/Documents/matrices';
fname = 'C:\Users\SiuLung\Downloads';

% enter case number
cn = 3;   % 1 = ctrl, 2 = dry, 3 = moist;

% enter grid spacing number
gs = 1;   % 1 = dx=125, 2 = dx=250, 5 = dx=500;

% get the corresponding filename
%---------------------------------------------
if cn == 1 && gs == 1     % ctrl case
filename = 'fmsfx_dx125_may08_4.mat';
filename1 = 'entm_thil_dx125_may08.mat';
filename2 = 'detm_thil_dx125_may08.mat';
end

if cn == 1 && gs == 2
filename = 'fmsfx_dx250_mar07_4.mat';
filename1 = 'entm_thil_dx250_mar07.mat';
filename2 = 'detm_thil_dx250_mar07.mat';
end

if cn == 1 && gs == 5
filename = 'fmsfx_dx500_feb27_4.mat';
filename1 = 'entm_thil_dx500_feb27.mat';
filename2 = 'detm_thil_dx500_feb27.mat';
end
%---------------------------------------------
if cn == 2 && gs == 1     % dry case
filename = 'fmsfx_dx125_dry_jul03_4.mat';
filename1 = 'entm_thil_dx125_dry_jul03.mat';
filename2 = 'detm_thil_dx125_dry_jul03.mat';
end

if cn == 2 && gs == 2
filename = 'fmsfx_dx250_dry_jul03_4.mat';
filename1 = 'entm_thil_dx250_dry_jul03.mat';
filename2 = 'detm_thil_dx250_dry_jul03.mat';
end

if cn == 2 && gs == 5
filename = 'fmsfx_dx500_dry_jul03_4.mat';
filename1 = 'entm_thil_dx500_dry_jul03.mat';
filename2 = 'detm_thil_dx500_dry_jul03.mat';
end
%---------------------------------------------
if cn == 3 && gs == 1     % moist case
filename = 'fmsfx_dx125_moist_jul03_4.mat';
filename1 = 'entm_thil_dx125_moist_jul03.mat';
filename2 = 'detm_thil_dx125_moist_jul03.mat';
end

if cn == 3 && gs == 2
filename = 'fmsfx_dx250_moist_jul03_4.mat';
filename1 = 'entm_thil_dx250_moist_jul03.mat';
filename2 = 'detm_thil_dx250_moist_jul03.mat';
end

if cn == 3 && gs == 5
filename = 'fmsfx_dx500_moist_jul03_4.mat';
filename1 = 'entm_thil_dx500_moist_jul03.mat';
filename2 = 'detm_thil_dx500_moist_jul03.mat';
end
%---------------------------------------------


fmsfx = load(fullfile(fname,filename));
entrain = load(fullfile(fname,filename1));

detm = { };

for seq = 6:6

ti = 0+(seq-1)*1; tf = 1+(seq-1)*1;

fieldname = sprintf('f%d%d',ti,tf);

fmf = fmsfx.fracmsfx.(fieldname).fmf;   % get frac. mass flux matrix
entm = entrain.entm_thil.(fieldname).ent;  % get frac. entrainment matrix (using theta_il)

% calculate fractional detrainment rate
det = entm - fmf;

% save work variable
detm_thil.(fieldname) = det;
save(fullfile(fname,filename2),'detm_thil');

end
