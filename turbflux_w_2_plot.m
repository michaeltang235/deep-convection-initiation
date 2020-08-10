close all
clear all

%----------------------------------------------------------------------------
% REMARKS:
% this code makes plots of filled contours of turb. flux of w (subgrid+resolved) and line
% contours of y-avg. w,
% it is used to see if updraft is widened in dx=125 and dx=500 when
% compared to dx=250 case
%----------------------------------------------------------------------------
tic 

% directory for dwdt_eqn.mat files
% fname = '/aos/home/mtang/Documents/matrices';
fname = 'C:\Users\SiuLung\Downloads';   % for dx = 250

% enter case number (1 = ctrl, 2 = dry, 3 = moist)
cn = 3;

% directory for output graphs
% fnameg = '/aos/home/mtang/Documents/graphs';
fnameg = 'C:\Users\SiuLung\Downloads';

% input range of z interested
zmax = 4;
zmin = 0;

% input plot number interested
% (1) - (5) are difference plots, (6) - (8) are full plots
% (1) = Th+Tv (subgrid flux from w-budget)
% (2) = Dw (subgrid flux from cm1 eqn. doc.)
% (3) = resolved turb. flux of w
% (4) = Th+Tv+resolved turb. flux of w
% (5) = Dw-resolved flux of w
% (6) = resolved turb. flux of w (full plots)
% (7) = Th+Tv-resolved turb. flux of w (fill plots)
% (8) = Dw-resolved flux of w (full plots)
% (9) = Profile of Dw-resolved flux of w (full plots)(avg. over. cen. 1 km)

% newly added option (difference plots)
% (10) = x-, y-, and z-comp. of resolved turb. flux of w on (3X3) subplots (n.a. in this version) 
% (20) = total (sub.+res.) turb. flux of w' in hori. and vert. on (3X2)subplots

pn = 20;

% get file names
if cn == 1
filename1 = 'turbflux_w_2_prcl_dx125_dec03.mat';
filename2 = 'turbflux_w_2_prcl_dx250_dec03.mat';
filename5 = 'turbflux_w_2_prcl_dx500_dec03.mat';
end
if cn == 2
filename1 = 'turbflux_w_2_dx125_dry_jul03.mat';
filename2 = 'turbflux_w_2_dx250_dry_jul03.mat';
filename5 = 'turbflux_w_2_dx500_dry_jul03.mat';
end
if cn == 3
filename1 = 'turbflux_w_2_dx125_moist_jul03.mat';
filename2 = 'turbflux_w_2_dx250_moist_jul03.mat';
filename5 = 'turbflux_w_2_dx500_moist_jul03.mat';
end

%-------------------------------------------------------------
% % get terms matrcies
terms1 = load(fullfile(fname,filename1));  % terms array  for dx=125
terms2 = load(fullfile(fname,filename2));  % terms array  for dx=250
terms5 = load(fullfile(fname,filename5));  % terms array  for dx=500

for seq = 6:6
%    
ti = 0+(seq-1)*1; tf = 1+(seq-1)*1;
% 
% for P = 1:3   % subplot position
% 
fieldname = sprintf('f%d%d',ti,tf);

% get xh and z matrices for all plots 
z = terms2.terms.(fieldname).z;  % in km
xh1 = terms1.terms.(fieldname).xh;   % in km for dx=125
xh2 = terms2.terms.(fieldname).xh;   % in km for dx=250
xh5 = terms5.terms.(fieldname).xh;   % in km for dx=500

% get xhii and xhfi from each case
xhii1 = terms1.terms.(fieldname).xhii;   % index for xmin for dx=125
xhii2 = terms2.terms.(fieldname).xhii;   % index for xmin for dx=250
xhii5 = terms5.terms.(fieldname).xhii;   % index for xmin for dx=500

xhfi1 = terms1.terms.(fieldname).xhfi;   % index for xmax for dx=125
xhfi2 = terms2.terms.(fieldname).xhfi;   % index for xmax for dx=250
xhfi5 = terms5.terms.(fieldname).xhfi;   % index for xmax for dx=500

% get index of zmax on z (half-level z)
zfi = find(abs(z-zmax)< 0.1/2);
zii = 1;

% get 2-d xh(xhii:xhfi) and z matrices for making contour plots
xh2d1n = repmat(xh1(xhii1:xhfi1),[1,zfi-zii+1]);  % dx=125
xh2d2n = repmat(xh2(xhii2:xhfi2),[1,zfi-zii+1]);  % dx=250
xh2d5n = repmat(xh5(xhii5:xhfi5),[1,zfi-zii+1]);  % dx=500

% add one lvl. at the bottom of xh2d (to remove white space as z starts at 0.5)
xh2d1 = [xh2d1n(:,1) xh2d1n];   % dx=125
xh2d2 = [xh2d2n(:,1) xh2d2n];   % dx=250
xh2d5 = [xh2d5n(:,1) xh2d5n];   % dx=500

z2d1n = repmat(z(zii:zfi)',[xhfi1-xhii1+1,1]);  % dx=125
z2d2n = repmat(z(zii:zfi)',[xhfi2-xhii2+1,1]);  % dx=250
z2d5n = repmat(z(zii:zfi)',[xhfi5-xhii5+1,1]);  % dx=500

% add one lvl. at the bottom of z (to remove white space as z starts at 0.5)

% create column matrix of zeros
z01 = zeros(size(z2d1n,1),1);   % dx=125
z02 = zeros(size(z2d2n,1),1);   % dx=250
z05 = zeros(size(z2d5n,1),1);   % dx=500

% concatenate z0 with z2d matrices
z2d1 = [z01 z2d1n];   % dx=125
z2d2 = [z02 z2d2n];   % dx=250
z2d5 = [z05 z2d5n];   % dx=500

if cn == 1
% get Th (hori. subgrid flux) of w-budget from arrays
hturbmrr1 = terms1.terms.(fieldname).hturbmrr;   % % time avg. of y-avg. of hori. turb. (subgrid flux from w-budget) for dx=125
hturbmrr2 = terms2.terms.(fieldname).hturbmrr;   % % time avg. of y-avg. of hori. turb. (subgrid flux from w-budget) for dx=250
hturbmrr5 = terms5.terms.(fieldname).hturbmrr;   % % time avg. of y-avg. of hori. turb. (subgrid flux from w-budget) for dx=500

% reduce hturbmrr matrices into chosen region
hturbmrr1rn = hturbmrr1(:,zii:zfi);   % for dx=125
hturbmrr2rn = hturbmrr2(:,zii:zfi);   % for dx=250
hturbmrr5rn = hturbmrr5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of hturbmrr2rn (to remove white space as z starts at 0.5)
hturbmrr1r = [hturbmrr1rn(:,1) hturbmrr1rn];   % for dx=125
hturbmrr2r = [hturbmrr2rn(:,1) hturbmrr2rn];   % for dx=250
hturbmrr5r = [hturbmrr5rn(:,1) hturbmrr5rn];   % for dx=500

% get Tv (vert. subgrid flux) of w-budget from arrays
vturbmrr1 = terms1.terms.(fieldname).vturbmrr;   % % time avg. of y-avg. of vert. turb. (subgrid flux from w-budget) for dx=125
vturbmrr2 = terms2.terms.(fieldname).vturbmrr;   % % time avg. of y-avg. of vert. turb. (subgrid flux from w-budget) for dx=250
vturbmrr5 = terms5.terms.(fieldname).vturbmrr;   % % time avg. of y-avg. of vert. turb. (subgrid flux from w-budget) for dx=500

% reduce vturbmrr matrices into chosen region
vturbmrr1rn = vturbmrr1(:,zii:zfi);   % for dx=125
vturbmrr2rn = vturbmrr2(:,zii:zfi);   % for dx=250
vturbmrr5rn = vturbmrr5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of vturbmrr2rn (to remove white space as z starts at 0.5)
vturbmrr1r = [vturbmrr1rn(:,1) vturbmrr1rn];   % for dx=125
vturbmrr2r = [vturbmrr2rn(:,1) vturbmrr2rn];   % for dx=250
vturbmrr5r = [vturbmrr5rn(:,1) vturbmrr5rn];   % for dx=500

% get Th+Tv (subgrid flux) of w-budget from arrays
thtvmrr1 = terms1.terms.(fieldname).thtvmrr;   % time avg. of y-avg. Th+Tv (subgrid flux) from w-budget for dx=125
thtvmrr2 = terms2.terms.(fieldname).thtvmrr;   % time avg. of y-avg. Th+Tv (subgrid flux) from w-budget for dx=250
thtvmrr5 = terms5.terms.(fieldname).thtvmrr;   % time avg. of y-avg. Th+Tv (subgrid flux) from w-budget for dx=500

% reduce thtvmrr matrices into chosen region
thtvmrr1rn = thtvmrr1(:,zii:zfi);   % for dx=125
thtvmrr2rn = thtvmrr2(:,zii:zfi);   % for dx=250
thtvmrr5rn = thtvmrr5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of thtvmrr2rn (to remove white space as z starts at 0.5)
thtvmrr1r = [thtvmrr1rn(:,1) thtvmrr1rn];   % for dx=125
thtvmrr2r = [thtvmrr2rn(:,1) thtvmrr2rn];   % for dx=250
thtvmrr5r = [thtvmrr5rn(:,1) thtvmrr5rn];   % for dx=500

% % get Dw (subgrid flux) from cm1 eqn. doc. from arrays
% sgfxm1 = terms1.terms.(fieldname).sgfxm;   % time avg. of y-avg. Dw (subgrid flux) from cm1 eqn. doc. for dx=125
% sgfxm2 = terms2.terms.(fieldname).sgfxm;   % time avg. of y-avg. Dw (subgrid flux) from cm1 eqn. doc. for dx=250
% sgfxm5 = terms5.terms.(fieldname).sgfxm;   % time avg. of y-avg. Dw (subgrid flux) from cm1 eqn. doc. for dx=500
% 
% % reduce sgfxm2 matrices into chosen region
% sgfxm1rn = sgfxm1(:,zii:zfi);   % for dx=125
% sgfxm2rn = sgfxm2(:,zii:zfi);   % for dx=250
% sgfxm5rn = sgfxm5(:,zii:zfi);   % for dx=500
% 
% % add one lvl. at the bottom of sgfxm2rn (to remove white space as z starts at 0.5)
% sgfxm1r = [sgfxm1rn(:,1) sgfxm1rn];   % for dx=125
% sgfxm2r = [sgfxm2rn(:,1) sgfxm2rn];   % for dx=250
% sgfxm5r = [sgfxm5rn(:,1) sgfxm5rn];   % for dx=500
end   % if cn == 1

% get resolved turb. flux of w from arrays
turbfxresm1 = terms1.terms.(fieldname).turbfxresm;   % time avg. of y-avg. resolved turb. flux of w for dx=125
turbfxresm2 = terms2.terms.(fieldname).turbfxresm;   % time avg. of y-avg. resolved turb. flux of w for dx=250
turbfxresm5 = terms5.terms.(fieldname).turbfxresm;   % time avg. of y-avg. resolved turb. flux of w for dx=500

% reduce turbfxm2 matrices into chosen region
turbfxresm1rn = turbfxresm1(:,zii:zfi);   % for dx=125
turbfxresm2rn = turbfxresm2(:,zii:zfi);   % for dx=250
turbfxresm5rn = turbfxresm5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of turbfxm2rn (to remove white space as z starts at 0.5)
turbfxresm1r = [turbfxresm1rn(:,1) turbfxresm1rn];   % for dx=125
turbfxresm2r = [turbfxresm2rn(:,1) turbfxresm2rn];   % for dx=250
turbfxresm5r = [turbfxresm5rn(:,1) turbfxresm5rn];   % for dx=500

% get Th+Tv-resolved flux of w (total flux) from arrays
thtvturbfxresm1 = terms1.terms.(fieldname).thtvturbfxresm;   % time avg. of y-avg. total turb. flux (Th+Tv-resolved flux) of w for dx=125
thtvturbfxresm2 = terms2.terms.(fieldname).thtvturbfxresm;   % time avg. of y-avg. total turb. flux (Th+Tv-resolved flux) of w for dx=250
thtvturbfxresm5 = terms5.terms.(fieldname).thtvturbfxresm;   % time avg. of y-avg. total turb. flux (Th+Tv-resolved flux) of w for dx=500

% reduce thtvturbfxm2 matrices into chosen region
thtvturbfxresm1rn = thtvturbfxresm1(:,zii:zfi);   % for dx=125
thtvturbfxresm2rn = thtvturbfxresm2(:,zii:zfi);   % for dx=250
thtvturbfxresm5rn = thtvturbfxresm5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of thtvturbfxm2rn (to remove white space as z starts at 0.5)
thtvturbfxresm1r = [thtvturbfxresm1rn(:,1) thtvturbfxresm1rn];   % for dx=125
thtvturbfxresm2r = [thtvturbfxresm2rn(:,1) thtvturbfxresm2rn];   % for dx=250
thtvturbfxresm5r = [thtvturbfxresm5rn(:,1) thtvturbfxresm5rn];   % for dx=500

% % get Dw-resolved flux of w (total flux) from arrays
% dwturbfxm1 = terms1.terms.(fieldname).dwturbfxm;   % time avg. of y-avg. total turb. flux (Dw-resolved flux) of w for dx=125
% dwturbfxm2 = terms2.terms.(fieldname).dwturbfxm;   % time avg. of y-avg. total turb. flux (Dw-resolved flux) of w for dx=250
% dwturbfxm5 = terms5.terms.(fieldname).dwturbfxm;   % time avg. of y-avg. total turb. flux (Dw-resolved flux) of w for dx=500
% 
% % reduce dwturbfxm2 matrices into chosen region
% dwturbfxm1rn = dwturbfxm1(:,zii:zfi);   % for dx=125
% dwturbfxm2rn = dwturbfxm2(:,zii:zfi);   % for dx=250
% dwturbfxm5rn = dwturbfxm5(:,zii:zfi);   % for dx=500
% 
% % add one lvl. at the bottom of dwturbfxm2rn (to remove white space as z starts at 0.5)
% dwturbfxm1r = [dwturbfxm1rn(:,1) dwturbfxm1rn];   % for dx=125
% dwturbfxm2r = [dwturbfxm2rn(:,1) dwturbfxm2rn];   % for dx=250
% dwturbfxm5r = [dwturbfxm5rn(:,1) dwturbfxm5rn];   % for dx=500

% get y-avg. w matrices from arrays
mywhrm1 = terms1.terms.(fieldname).mywhrm;   % time avg. of y-avg. w for dx=125
mywhrm2 = terms2.terms.(fieldname).mywhrm;   % time avg. of y-avg. w for dx=250
mywhrm5 = terms5.terms.(fieldname).mywhrm;   % time avg. of y-avg. w for dx=500

% reduce mywhrm matrices into chosen region
mywhrm1rn = mywhrm1(:,zii:zfi);   % for dx=125
mywhrm2rn = mywhrm2(:,zii:zfi);   % for dx=250
mywhrm5rn = mywhrm5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of mywhrm1rn (to remove white space as z starts at 0.5)
mywhrm1r = [mywhrm1rn(:,1) mywhrm1rn];   % for dx=125
mywhrm2r = [mywhrm2rn(:,1) mywhrm2rn];   % for dx=250
mywhrm5r = [mywhrm5rn(:,1) mywhrm5rn];   % for dx=500

%-------------------
% get hori. comp. of resolved flux of w (total flux) from arrays
turbfxreshm1 = terms1.terms.(fieldname).turbfxreshm;   % time avg. of y-avg. of hori.-compo. of resolved turb. flux of w for dx=125
turbfxreshm2 = terms2.terms.(fieldname).turbfxreshm;   % time avg. of y-avg. of hori.-compo. of resolved turb. flux of w for dx=250
turbfxreshm5 = terms5.terms.(fieldname).turbfxreshm;   % time avg. of y-avg. of hori.-compo. of resolved turb. flux of w for dx=500

% reduce turbfxreshm2 matrices into chosen region
turbfxreshm1rn = turbfxreshm1(:,zii:zfi);   % for dx=125
turbfxreshm2rn = turbfxreshm2(:,zii:zfi);   % for dx=250
turbfxreshm5rn = turbfxreshm5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of turbfxreshm2rn (to remove white space as z starts at 0.5)
turbfxreshm1r = [turbfxreshm1rn(:,1) turbfxreshm1rn];   % for dx=125
turbfxreshm2r = [turbfxreshm2rn(:,1) turbfxreshm2rn];   % for dx=250
turbfxreshm5r = [turbfxreshm5rn(:,1) turbfxreshm5rn];   % for dx=500
%-------------------
% get vert.-comp. of resolved flux of w (total flux) from arrays
turbfxresvm1 = terms1.terms.(fieldname).turbfxresvm;   % time avg. of y-avg. of vert.-compo. of resolved turb. flux of w for dx=125
turbfxresvm2 = terms2.terms.(fieldname).turbfxresvm;   % time avg. of y-avg. of vert.-compo. of resolved turb. flux of w for dx=250
turbfxresvm5 = terms5.terms.(fieldname).turbfxresvm;   % time avg. of y-avg. of vert.-compo. of resolved turb. flux of w for dx=500

% reduce turbfxresvm2 matrices into chosen region
turbfxresvm1rn = turbfxresvm1(:,zii:zfi);   % for dx=125
turbfxresvm2rn = turbfxresvm2(:,zii:zfi);   % for dx=250
turbfxresvm5rn = turbfxresvm5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of turbfxresvm2rn (to remove white space as z starts at 0.5)
turbfxresvm1r = [turbfxresvm1rn(:,1) turbfxresvm1rn];   % for dx=125
turbfxresvm2r = [turbfxresvm2rn(:,1) turbfxresvm2rn];   % for dx=250
turbfxresvm5r = [turbfxresvm5rn(:,1) turbfxresvm5rn];   % for dx=500

%--------------------------------------------------------------------------------
% set axes font size
axfs = 12;   % axes font size
tfs = 14;   % title font size
tefs = 11;   % text font size (case number, e.g., dx=250 m)
lbfs = 14;   % panel label font size
cbfs = 11;   % colorbar font size

xmina = -3;   % xmin for all plots
xmaxa = 3;   % xmax for all plots
xticka = 1.5;   % xticks for all plots

pbra = 1/3;   % apsect ratio for plots, y:x ratio

clv = -0.006:0.006/10:0.006;   % for total turb. flux of w contours

% if pn == 1
%     clv = -7.5e-4:1e-4:7.5e-4;
%     dclv = -3.5e-4:1e-4:3.5e-4;
% end

if pn == 1
    clv = -0.0014:0.0004:0.0014;
    dclv = -6.5e-4:1e-4:6.5e-4;
end

if pn == 2
    clv = -0.0014:0.0004:0.0014;
    dclv = -6.5e-4:1e-4:6.5e-4;
end

if pn == 3
    clv = -0.0055:0.001:0.0055;
    dclv = -0.00225:0.001/2:0.00225;
end

if pn == 4
    clv = -0.0055:0.0002:0.0055;
    dclv = 0.5.*clv;
%     clv = -0.0055:0.001:0.0055;
%     dclv = -0.00225:0.001/2:0.00225;
end

if pn == 5
    clv = -0.0055:0.001:0.0055;
    dclv = -0.00225:0.001/2:0.00225;
end

if pn == 6 || pn == 7 || pn == 8
%     clv = -0.0065:0.001:0.0065;
    clv = -0.0066:0.0004:0.0066;
end

if pn == 10 || pn == 11 || pn == 12 || pn == 20
    clv = -0.0055:0.0002:0.0055;
    dclv = 0.5.*clv;
%     clv = -0.0055:0.001:0.0055;
%     dclv = -0.00225:0.001/2:0.00225;
end

% 
% llvp = 0:0.5:4;   % for positive y-avg. w contours

llvp = 0.25:0.5:3.75;   % for positive y-avg. w contours
llvn = -1.75:0.5:-0.25;   % for negative y-avg. w contours   % for negative y-avg. w contours

dllvp = 0.1:0.2:0.5;   % for diff. in positive y-avg. w contours
dllvn = -0.9:0.2:-0.1;   % for diff. in negative y-avg. w contours

%--------------------------------------------------------------------------------
% OPEN FIGURE WINDOW

f = figure('units','normalized','outerposition',[0 0 1 1]);
% f = figure('units','normalized');
%---------------------------------------------------------

% (I) MAKE PLOT OF Th+Tv (SUBGRID FLUX OF W FROM W-BUDGET)
if pn == 1
    
% FOR dx=125: GET DIFFEENCE PLOT

% interpolate thtvmrr1r and myw from dx=125 to dx=250 case
thtvmrr1r2 = 0.5.*(thtvmrr1r(1:2:end-1,:)+ thtvmrr1r(2:2:end,:));
mywhrm1r2 = 0.5.*(mywhrm1r(1:2:end-1,:) + mywhrm1r(2:2:end,:));

% get difference btw. interpolated thtvmrr1r2 from thtvmrr2r (mywhrm1r2
% from mywhrm2r)
diffthtvmrr12 = thtvmrr1r2 - thtvmrr2r;
diffmywhrm12 = mywhrm1r2 - mywhrm2r;

subplot(3,1,1)

[C11, h11] = contourf(xh2d2,z2d2,diffthtvmrr12,dclv); hold on 
[C21, h21] = contour(xh2d2,z2d2,diffmywhrm12,dllvp); hold on
[C31, h31] = contour(xh2d2,z2d2,diffmywhrm12,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h11.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
% cb.TickLabels = dclv;
% cb.Ticks = dclv;
cb.Ticks = -0.0006:0.0002:0.0006;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title([ '\(\overline{\textrm{T}_h + \textrm{T}_v}\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1 = gca;   % for dx=125 m
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(2.5,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(0.85,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% END for dx=125
%---------------------------------------------------------
% FOR dx=250 

subplot(3,1,2)

[C12, h12] = contourf(xh2d2,z2d2,thtvmrr2r,clv); hold on 
[C22, h22] = contour(xh2d2,z2d2,mywhrm2r,llvp); hold on
[C32, h32] = contour(xh2d2,z2d2,mywhrm2r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h12.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
cb.Ticks = clv;
% cb.Ruler.Exponent = -3;
% cb.Ticks = clv([1 3 5 7 10 12 14 16]);
% cb.Ticks = -0.0007:0.0002:0.0007;
% cb.Ticks = clv(1:2:end);

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2 = gca;   % for dx=250 m
ax2.XLim = [xmin,xmax];
ax2.XTick = [xmin:xtick:xmax];

ax2.YLim = [ymin,ymax];
ax2.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(2,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=250

%---------------------------------------------------------

% FOR dx = 500 m: GET DIFFERENCE PLOT

% interpolate thtvmrr2r and myw from dx=250 to dx=500 case
thtvmrr2r5 = 0.5.*(thtvmrr2r(1:2:end-1,:)+ thtvmrr2r(2:2:end,:));
mywhrm2r5 = 0.5.*(mywhrm2r(1:2:end-1,:) + mywhrm2r(2:2:end,:));

% get difference btw. interpolated tturbfxwm5r from tturbfxwm2r5 (mywhrm5r
% from mywhrm2r5)
diffthtvmrr52 = thtvmrr5r - thtvmrr2r5;
diffmywhrm52 = mywhrm5r - mywhrm2r5;

% To prevent white patches formed on contour plot because of entries in
% diffthtvmrr52 exceeding limits set by dclv,
% find those entries and set them as the limit of dclv

% initialize new diffthtvmrr52
diffthtvmrr52new = diffthtvmrr52;

if find(diffthtvmrr52 > dclv(end)) ~= 0   % for upper limit of dclv
    diffthtvmrr52new(find(diffthtvmrr52 > dclv(end))) = dclv(end);
end

if find(diffthtvmrr52 < dclv(1)) ~= 0   % for lower limit of dclv
    diffthtvmrr52new(find(diffthtvmrr52 < dclv(1))) = dclv(1);
end

subplot(3,1,3)

[C15, h15] = contourf(xh2d5,z2d5,diffthtvmrr52new,dclv); hold on 
[C25, h25] = contour(xh2d5,z2d5,diffmywhrm52,dllvp); hold on
[C35, h35] = contour(xh2d5,z2d5,diffmywhrm52,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h15.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
% cb.Ticks = dclv;
cb.Ticks = -0.0006:0.0002:0.0006;

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5 = gca;   % for dx=500 m
ax5.XLim = [xmin,xmax];
ax5.XTick = [xmin:xtick:xmax];

ax5.YLim = [ymin,ymax];
ax5.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(0.85,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=500
%---------------------------------------------------------
end   % end if pn == 1   

% END MAKING PLOT OF (Th+Tv, subgrid flux of w from w-budget)
%----------------------------------------------------------------------------

% (II) MAKE PLOT OF Dw (SUBGRID FLUX FROM CM1 EQN. DOC.)
if pn == 2

% FOR dx=125: GET DIFFEENCE PLOT

% interpolate sgfxm1r and myw from dx=125 to dx=250 case
sgfxm1r2 = 0.5.*(sgfxm1r(1:2:end-1,:)+ sgfxm1r(2:2:end,:));
mywhrm1r2 = 0.5.*(mywhrm1r(1:2:end-1,:) + mywhrm1r(2:2:end,:));

% get difference btw. interpolated sgfxm1r2 from sgfxm2r (mywhrm1r2
% from mywhrm2r)
diffsgfxm12 = sgfxm1r2 - sgfxm2r;
diffmywhrm12 = mywhrm1r2 - mywhrm2r;

subplot(3,1,1)

[C11, h11] = contourf(xh2d2,z2d2,diffsgfxm12,dclv); hold on 
[C21, h21] = contour(xh2d2,z2d2,diffmywhrm12,dllvp); hold on
[C31, h31] = contour(xh2d2,z2d2,diffmywhrm12,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h11.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
% cb.Ticks = dclv;
% cb.Ticks = dclv([1 3 5 7 8 10 12 14]);
cb.Ticks = -0.0006:0.0002:0.0006;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title([ '\(\overline{\textrm{D}_w}\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1 = gca;   % for dx=125 m
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(2.5,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(0.85,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% END for dx=125
%---------------------------------------------------------
% FOR dx=250 

subplot(3,1,2)

[C12, h12] = contourf(xh2d2,z2d2,sgfxm2r,clv); hold on 
[C22, h22] = contour(xh2d2,z2d2,mywhrm2r,llvp); hold on
[C32, h32] = contour(xh2d2,z2d2,mywhrm2r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h12.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
cb.Ticks = clv;
% cb.Ruler.Exponent = -3;
% cb.Ticks = clv([1 3 5 7 10 12 14 16]);
% cb.Ticks = clv(1:4:end);

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2 = gca;   % for dx=250 m
ax2.XLim = [xmin,xmax];
ax2.XTick = [xmin:xtick:xmax];

ax2.YLim = [ymin,ymax];
ax2.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(2,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=250

%---------------------------------------------------------

% FOR dx = 500 m: GET DIFFERENCE PLOT

% interpolate sgfxm2r and myw from dx=250 to dx=500 case
sgfxm2r5 = 0.5.*(sgfxm2r(1:2:end-1,:)+ sgfxm2r(2:2:end,:));
mywhrm2r5 = 0.5.*(mywhrm2r(1:2:end-1,:) + mywhrm2r(2:2:end,:));

% get difference btw. interpolated sgfxm5r from sgfxm2r5 (mywhrm5r
% from mywhrm2r5)
diffsgfxm52 = sgfxm5r - sgfxm2r5;
diffmywhrm52 = mywhrm5r - mywhrm2r5;

% To prevent white patches formed on contour plot because of entries in
% diffsgfxm52 exceeding limits set by dclv,
% find those entries and set them as the limit of dclv

% initialize new diffsgfxm52
diffsgfxm52new = diffsgfxm52;

if find(diffsgfxm52 > dclv(end)) ~= 0   % for upper limit of dclv
    diffsgfxm52new(find(diffsgfxm52 > dclv(end))) = dclv(end);
end

if find(diffsgfxm52 < dclv(1)) ~= 0   % for lower limit of dclv
    diffsgfxm52new(find(diffsgfxm52 < dclv(1))) = dclv(1);
end

subplot(3,1,3)

[C15, h15] = contourf(xh2d5,z2d5,diffsgfxm52new,dclv); hold on 
[C25, h25] = contour(xh2d5,z2d5,diffmywhrm52,dllvp); hold on
[C35, h35] = contour(xh2d5,z2d5,diffmywhrm52,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h15.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
% cb.Ticks = dclv;
% cb.Ticks = dclv([1 3 5 7 8 10 12 14]);
cb.Ticks = -0.0006:0.0002:0.0006;

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5 = gca;   % for dx=500 m
ax5.XLim = [xmin,xmax];
ax5.XTick = [xmin:xtick:xmax];

ax5.YLim = [ymin,ymax];
ax5.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(0.85,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=500
%---------------------------------------------------------
end   % end if pn == 2 

% END MAKING PLOT OF (Dw, subgrid flux of w from cm1 eqn. doc.)
%----------------------------------------------------------------------------

% (III) MAKE PLOT OF RESOLVED TURB. FLUX OF W
if pn == 3

% FOR dx=125: GET DIFFEENCE PLOT

% interpolate turbfxm1r and myw from dx=125 to dx=250 case
turbfxm1r2 = 0.5.*(turbfxm1r(1:2:end-1,:)+ turbfxm1r(2:2:end,:));
mywhrm1r2 = 0.5.*(mywhrm1r(1:2:end-1,:) + mywhrm1r(2:2:end,:));

% get difference btw. interpolated turbfxm1r2 from turbfxm2r (mywhrm1r2
% from mywhrm2r)
diffturbfxm12 = -turbfxm1r2 - (-turbfxm2r);   % - turb. flux|125 - (-turb. flux|250)
diffmywhrm12 = mywhrm1r2 - mywhrm2r;

subplot(3,1,1)

[C11, h11] = contourf(xh2d2,z2d2,diffturbfxm12,dclv); hold on 
[C21, h21] = contour(xh2d2,z2d2,diffmywhrm12,dllvp); hold on
[C31, h31] = contour(xh2d2,z2d2,diffmywhrm12,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h11.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title([ '\(-\nabla\cdot (\overline{\textrm{\textbf{u}}^{\prime} w^{\prime}})\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% title([ '\(\overline{\textrm{resolved turb. flux}}\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1 = gca;   % for dx=125 m
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(2.5,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(0.85,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% END for dx=125
%---------------------------------------------------------
% FOR dx=250 

subplot(3,1,2)

[C12, h12] = contourf(xh2d2,z2d2,-turbfxm2r,clv); hold on 
[C22, h22] = contour(xh2d2,z2d2,mywhrm2r,llvp); hold on
[C32, h32] = contour(xh2d2,z2d2,mywhrm2r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h12.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
% cb.Ruler.Exponent = -3;
% cb.Ticks = clv([1 3 5 8 10 12]);
cb.Ticks = -0.005:0.0025:0.005;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2 = gca;   % for dx=250 m
ax2.XLim = [xmin,xmax];
ax2.XTick = [xmin:xtick:xmax];

ax2.YLim = [ymin,ymax];
ax2.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(2,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=250

%---------------------------------------------------------

% FOR dx = 500 m: GET DIFFERENCE PLOT

% interpolate turbfxm2r and myw from dx=250 to dx=500 case
dwturbfxm2r5 = 0.5.*(turbfxm2r(1:2:end-1,:)+ turbfxm2r(2:2:end,:));
mywhrm2r5 = 0.5.*(mywhrm2r(1:2:end-1,:) + mywhrm2r(2:2:end,:));

% get difference btw. interpolated turbfxm5r from turbfxm2r5 (mywhrm5r
% from mywhrm2r5)
diffturbfxm52 = -turbfxm5r - (-dwturbfxm2r5);   % -turb. flux|500 - (-turb. flux|250)
diffmywhrm52 = mywhrm5r - mywhrm2r5;

% To prevent white patches formed on contour plot because of entries in
% diffturbfxm52 exceeding limits set by dclv,
% find those entries and set them as the limit of dclv

% initialize new diffturbfxm52
diffturbfxm52new = diffturbfxm52;

if find(diffturbfxm52 > dclv(end)) ~= 0   % for upper limit of dclv
    diffturbfxm52new(find(diffturbfxm52 > dclv(end))) = dclv(end);
end

if find(diffturbfxm52 < dclv(1)) ~= 0   % for lower limit of dclv
    diffturbfxm52new(find(diffturbfxm52 < dclv(1))) = dclv(1);
end

subplot(3,1,3)

[C15, h15] = contourf(xh2d5,z2d5,diffturbfxm52new,dclv); hold on 
[C25, h25] = contour(xh2d5,z2d5,diffmywhrm52,dllvp); hold on
[C35, h35] = contour(xh2d5,z2d5,diffmywhrm52,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h15.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5 = gca;   % for dx=500 m
ax5.XLim = [xmin,xmax];
ax5.XTick = [xmin:xtick:xmax];

ax5.YLim = [ymin,ymax];
ax5.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(0.85,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=500
%---------------------------------------------------------
end   % end if pn == 3

% END MAKING PLOT OF (resolved turb. flux of w
%----------------------------------------------------------------------------

% (IV) MAKE PLOT OF Th+Tv-RESOLVED TURB. FLUX OF W
if pn == 4

% FOR dx=125: GET DIFFEENCE PLOT

% % interpolate thtvturbfxm1r and myw from dx=125 to dx=250 case
% thtvturbfxresm1r2 = 0.5.*(thtvturbfxresm1r(1:2:end-1,:)+ thtvturbfxresm1r(2:2:end,:));
% mywhrm1r2 = 0.5.*(mywhrm1r(1:2:end-1,:) + mywhrm1r(2:2:end,:));
% 
% % get difference btw. interpolated thtvturbfxm1r2 from thtvturbfxm2r (mywhrm1r2
% % from mywhrm2r)
% diffthtvturbfxresm12 = thtvturbfxresm1r2 - thtvturbfxresm2r;
% diffmywhrm12 = mywhrm1r2 - mywhrm2r;

%-----------------------
% use interp1 to interpolate grid points from 250-m to 125-m case
thtvturbfxresm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));   % y-avg.total (subgrid+res.) turb. flux of w' from dx=250 to dx=125
mywhrm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));  % y-avg. w from dx=250 to dx=125

for k = 1:size(mywhrm1r,2)
thtvturbfxresm2r1(:,k) = interp1(xh2(xhii2:xhfi2),thtvturbfxresm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
mywhrm2r1(:,k) = interp1(xh2(xhii2:xhfi2),mywhrm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
end

% get diff. between 125-m and 250-m cases
diffthtvturbfxresm12 = thtvturbfxresm1r - thtvturbfxresm2r1;
diffmywhrm12 = mywhrm1r - mywhrm2r1;   % diff. in w (dx=500-dx=250)

%-----------------------

subplot(3,1,1)

[C11, h11] = contourf(xh2d1,z2d1,diffthtvturbfxresm12,dclv); hold on 
[C21, h21] = contour(xh2d1,z2d1,diffmywhrm12,dllvp); hold on
[C31, h31] = contour(xh2d1,z2d1,diffmywhrm12,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h11.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb1 = colorbar;
cb1.Position(1) = 0.9;
cb1.FontSize = cbfs;
cb1.TickLabelInterpreter = 'latex';

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}\big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

if cn == 1
title('\(\overline{\textrm{T}}_{h,\mathrm{sgs}} + \overline{\textrm{T}}_{v,\mathrm{sgs}} - \frac{1}{\rho_{0}}\nabla\cdot(\rho_{0}\overline{\textrm{\textbf{u}}^{\prime}}\overline{w^{\prime}})\)','Interpreter','latex','FontSize',tfs);
end
if cn ~= 1
title('\(- \frac{1}{\rho_{0}}\nabla\cdot(\rho_{0}\overline{\textrm{\textbf{u}}^{\prime}}\overline{w^{\prime}})\)','Interpreter','latex','FontSize',tfs);
end
% title('\(\overline{\textrm{T}_{h,\mathrm{sgs}} + \textrm{T}_{v,\mathrm{sgs}} - \frac{1}{\rho_{0}}\nabla\cdot(\rho_{0}\textrm{\textbf{u}}^{\prime}w^{\prime})}\)','Interpreter','latex','FontSize',tfs);

% title(['\(\textrm{T}_h + \textrm{T}_v - \textrm{resolved turb. flux} \), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1 = gca;   % for dx=125 m
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];
ax1.FontSize = axfs;
ax1.Title.FontSize = tfs;
ax1.TickLabelInterpreter = 'latex';

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(-2.90,2.6,'(a)', 'Interpreter', 'latex','FontSize',lbfs,'FontWeight','bold');   % upper left corner

text(0.70,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% END for dx=125
%---------------------------------------------------------
% FOR dx=250 

subplot(3,1,2)
[C12, h12] = contourf(xh2d2,z2d2,thtvturbfxresm2r,clv); hold on 
[C22, h22] = contour(xh2d2,z2d2,mywhrm2r,llvp); hold on
[C32, h32] = contour(xh2d2,z2d2,mywhrm2r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h12.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb2 = colorbar;
cb2.Ticks = -0.005:0.0025:0.005;
cb2.Position(1) = 0.9;
cb2.FontSize = cbfs;
cb2.TickLabelInterpreter = 'latex';
% cb.Ruler.Exponent = -3;
% cb.Ticks = clv(1:2:end);

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z} \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2 = gca;   % for dx=250 m
ax2.XLim = [xmin,xmax];
ax2.XTick = [xmin:xtick:xmax];

ax2.YLim = [ymin,ymax];
ax2.YTick = [ymin:ytick:ymax];
ax2.FontSize = axfs;
ax2.TickLabelInterpreter = 'latex';


pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(-2.90,2.6,'(b)', 'Interpreter', 'latex','FontSize',lbfs,'FontWeight','bold');   % upper left corner

text(1.9,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

if cn == 1 
    text(2.35,0.25,'CTRL', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
end
if cn == 2
    text(2.45,0.25,'DRY', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
end
if cn == 3 
    text(2.27,0.25,'MOIST', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
end

% end for dx=250

%---------------------------------------------------------

% FOR dx = 500 m: GET DIFFERENCE PLOT

% % interpolate thtvturbfxm5r and myw from dx=250 to dx=500 case
% thtvturbfxresm2r5 = 0.5.*(thtvturbfxresm2r(1:2:end-1,:)+ thtvturbfxresm2r(2:2:end,:));
% mywhrm2r5 = 0.5.*(mywhrm2r(1:2:end-1,:) + mywhrm2r(2:2:end,:));
% 
% % get difference btw. interpolated thtvturbfxm5r from thtvturbfxm2r5 (mywhrm5r
% % from mywhrm2r5)
% diffthtvturbfxresm52 = thtvturbfxresm5r - thtvturbfxresm2r5;
% diffmywhrm52 = mywhrm5r - mywhrm2r5;


%-----------------------
% use interp1 to interpolate grid points from 500-m to 250-m case
thtvturbfxresm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));   % y-avg.total (subgrid+res.) turb. flux of w' from dx=500 to dx=250
mywhrm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm2r,2)
thtvturbfxresm5r2(:,k) = interp1(xh5(xhii5:xhfi5),thtvturbfxresm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
mywhrm5r2(:,k) = interp1(xh5(xhii5:xhfi5),mywhrm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
end

% get diff. between 500-m and 250-m cases
diffthtvturbfxresm52 = thtvturbfxresm5r2 - thtvturbfxresm2r;
diffmywhrm52 = mywhrm5r2 - mywhrm2r;   % diff. in w (dx=500-dx=250)

%-----------------------


% To prevent white patches formed on contour plot because of entries in
% diffthtvturbfxm52 exceeding limits set by dclv,
% find those entries and set them as the limit of dclv

% initialize new diffthtvturbfxm52
diffthtvturbfxresm52new = diffthtvturbfxresm52;

if find(diffthtvturbfxresm52 > dclv(end)) ~= 0   % for upper limit of dclv
    diffthtvturbfxresm52new(find(diffthtvturbfxresm52 > dclv(end))) = dclv(end);
end

if find(diffthtvturbfxresm52 < dclv(1)) ~= 0   % for lower limit of dclv
    diffthtvturbfxresm52new(find(diffthtvturbfxresm52 < dclv(1))) = dclv(1);
end

subplot(3,1,3)

[C15, h15] = contourf(xh2d2,z2d2,diffthtvturbfxresm52new,dclv); hold on 
[C25, h25] = contour(xh2d2,z2d2,diffmywhrm52,dllvp); hold on
[C35, h35] = contour(xh2d2,z2d2,diffmywhrm52,dllvn); hold off

% [C225, h225] = contour(xh2d5,z2d5,mywhrm2r5,llvp); hold on
% [C335, h335] = contour(xh2d5,z2d5,mywhrm2r5,llvn); hold on
% 
% h225.LineColor = 'b';   % set contour lines at black

% [C25, h25] = contour(xh2d5,z2d5,mywhrm5r,llvp); hold on
% [C35, h35] = contour(xh2d5,z2d5,mywhrm5r,llvn); hold on

grid on 
box on

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h15.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb5 = colorbar;
cb5.Position(1) = 0.9;
cb5.FontSize = cbfs;
cb5.TickLabelInterpreter = 'latex';
% cb.Ticks = -0.002:0.001:0.002;

xlabel('\textit{x} \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z} \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5 = gca;   % for dx=500 m
ax5.XLim = [xmin,xmax];
ax5.XTick = [xmin:xtick:xmax];

ax5.YLim = [ymin,ymax];
ax5.YTick = [ymin:ytick:ymax];
ax5.FontSize = axfs;
ax5.TickLabelInterpreter = 'latex';


pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(-2.90,2.6,'(c)', 'Interpreter', 'latex','FontSize',lbfs,'FontWeight','bold');   % upper left corner

text(0.70,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=500
%---------------------------------------------------------

% SET SUBPLOT POSITIONS

% reduce vert. spacing between subplots
ax2.Position(2) = ax5.Position(2) + ax5.Position(4) + 0.0425;
ax1.Position(2) = ax2.Position(2) + ax2.Position(4) + 0.0425;

% set colorbar width
cb1.Position(3) = 0.5*cb1.Position(3);
cb2.Position(3) = 0.5*cb2.Position(3);
cb5.Position(3) = 0.5*cb5.Position(3);

% reduce hori. spacing between subplots and colorbars
cb1.Position(1) = ax1.Position(1) + ax1.Position(3) - 0.2175;
cb2.Position(1) = ax2.Position(1) + ax2.Position(3) - 0.2175;
cb5.Position(1) = ax5.Position(1) + ax5.Position(3) - 0.2175;

% make colorbars level as subplots
cb1.Position(2) = ax1.Position(2);
cb2.Position(2) = ax2.Position(2);
cb5.Position(2) = ax5.Position(2);

% DONE SETTING SUBPLOT POSITION
%---------------------------------------------------------

end   % end if pn == 4

% END MAKING PLOT OF (Th+Tv+resolved turb. flux of w)
%----------------------------------------------------------------------------

% (V) MAKE PLOT OF Dw-RESOLVED TURB. FLUX OF W
if pn == 5

% FOR dx=125: GET DIFFEENCE PLOT

% interpolate dwturbfxm1r and myw from dx=125 to dx=250 case
dwturbfxm1r2 = 0.5.*(dwturbfxm1r(1:2:end-1,:)+ dwturbfxm1r(2:2:end,:));
mywhrm1r2 = 0.5.*(mywhrm1r(1:2:end-1,:) + mywhrm1r(2:2:end,:));

% get difference btw. interpolated dwturbfxm1r2 from dwturbfxm2r (mywhrm1r2
% from mywhrm2r)
diffdwturbfxm12 = dwturbfxm1r2 - dwturbfxm2r;
diffmywhrm12 = mywhrm1r2 - mywhrm2r;

subplot(3,1,1)

[C11, h11] = contourf(xh2d2,z2d2,diffdwturbfxm12,dclv); hold on 
[C21, h21] = contour(xh2d2,z2d2,diffmywhrm12,dllvp); hold on
[C31, h31] = contour(xh2d2,z2d2,diffmywhrm12,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h11.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\overline{\textrm{D}_w - \nabla\cdot (\overline{\textrm{\textbf{u}}^{\prime} w^{\prime}})}\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% title([ '\(\overline{\textrm{D}_w - \textrm{resolved turb. flux}}\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1 = gca;   % for dx=125 m
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(2.5,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(0.85,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% END for dx=125
%---------------------------------------------------------
% FOR dx=250 

subplot(3,1,2)
[C12, h12] = contourf(xh2d2,z2d2,dwturbfxm2r,clv); hold on 
[C22, h22] = contour(xh2d2,z2d2,mywhrm2r,llvp); hold on
[C32, h32] = contour(xh2d2,z2d2,mywhrm2r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h12.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
% cb.Ruler.Exponent = -3;
cb.Ticks = -0.005:0.0025:0.005;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2 = gca;   % for dx=250 m
ax2.XLim = [xmin,xmax];
ax2.XTick = [xmin:xtick:xmax];

ax2.YLim = [ymin,ymax];
ax2.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(2,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=250

%---------------------------------------------------------

% FOR dx = 500 m: GET DIFFERENCE PLOT

% interpolate dwturbfxm2r and myw from dx=250 to dx=500 case
dwturbfxm2r5 = 0.5.*(dwturbfxm2r(1:2:end-1,:)+ dwturbfxm2r(2:2:end,:));
mywhrm2r5 = 0.5.*(mywhrm2r(1:2:end-1,:) + mywhrm2r(2:2:end,:));

% get difference btw. interpolated dwturbfxm5r from dwturbfxm2r5 (mywhrm5r
% from mywhrm2r5)
diffdwturbfxm52 = dwturbfxm5r - dwturbfxm2r5;
diffmywhrm52 = mywhrm5r - mywhrm2r5;

% To prevent white patches formed on contour plot because of entries in
% diffdwturbfxm52 exceeding limits set by dclv,
% find those entries and set them as the limit of dclv

% initialize new diffdwturbfxm52
diffdwturbfxm52new = diffdwturbfxm52;

if find(diffdwturbfxm52 > dclv(end)) ~= 0   % for upper limit of dclv
    diffdwturbfxm52new(find(diffdwturbfxm52 > dclv(end))) = dclv(end);
end

if find(diffdwturbfxm52 < dclv(1)) ~= 0   % for lower limit of dclv
    diffdwturbfxm52new(find(diffdwturbfxm52 < dclv(1))) = dclv(1);
end

subplot(3,1,3)

[C15, h15] = contourf(xh2d5,z2d5,diffdwturbfxm52new,dclv); hold on 
[C25, h25] = contour(xh2d5,z2d5,diffmywhrm52,dllvp); hold on
[C35, h35] = contour(xh2d5,z2d5,diffmywhrm52,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h15.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5 = gca;   % for dx=500 m
ax5.XLim = [xmin,xmax];
ax5.XTick = [xmin:xtick:xmax];

ax5.YLim = [ymin,ymax];
ax5.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(0.85,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=500
%---------------------------------------------------------
end   % end if pn == 5

% END MAKING PLOT OF (Dw-resolved turb. flux of w)
%----------------------------------------------------------------------------

% (VI) MAKE PLOT OF RESOLVED TURB. FLUX OF W (FULL PLOTS)
if pn == 6

% FOR dx=125: GET FULL PLOT (NOT DIFFERENCE PLOT)

subplot(3,1,1)

[C11, h11] = contourf(xh2d1,z2d1,-turbfxm1r,clv); hold on 
[C21, h21] = contour(xh2d1,z2d1,mywhrm1r,llvp); hold on
[C31, h31] = contour(xh2d1,z2d1,mywhrm1r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h11.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
cb.Ticks = -0.006:0.003:0.006;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title([ '\(-\nabla\cdot (\overline{\textrm{\textbf{u}}^{\prime} w^{\prime}})\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% title([ '\(\overline{\textrm{resolved turb. flux}}\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1 = gca;   % for dx=125 m
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(2.5,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% text(0.85,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(2,2.7,'\(\Delta=125\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% END for dx=125
%---------------------------------------------------------
% FOR dx=250 

subplot(3,1,2)

[C12, h12] = contourf(xh2d2,z2d2,-turbfxm2r,clv); hold on 
[C22, h22] = contour(xh2d2,z2d2,mywhrm2r,llvp); hold on
[C32, h32] = contour(xh2d2,z2d2,mywhrm2r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h12.LineColor = 'none';   % set dudx contour lines as colorless

% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap

blue=[0,0,1];
red=[1,0,0];
white=[1,1,1];
ncolors=15;
nc2=0.5*(ncolors-1);
for nc=1:ncolors
  if (nc <= nc2+1)
    bluered(nc,:)=(nc-1)/nc2*white+(nc2-nc+1)/nc2*blue;
  else
    bluered(nc,:)=(nc-nc2-1)/nc2*red+(ncolors-nc)/nc2*white;
  end
end

colormap(bluered);   %insert chosen colormap


cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
% cb.Ruler.Exponent = -3;
% cb.Ticks = clv([1 3 5 8 10 12]);
cb.Ticks = -0.006:0.003:0.006;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2 = gca;   % for dx=250 m
ax2.XLim = [xmin,xmax];
ax2.XTick = [xmin:xtick:xmax];

ax2.YLim = [ymin,ymax];
ax2.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(2,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=250

%---------------------------------------------------------

% FOR dx = 500 m: GET FULL PLOT (NOT DIFFERENCE PLOT)

subplot(3,1,3)

[C15, h15] = contourf(xh2d5,z2d5,-turbfxm5r,clv); hold on 
[C25, h25] = contour(xh2d5,z2d5,mywhrm5r,llvp); hold on
[C35, h35] = contour(xh2d5,z2d5,mywhrm5r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h15.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
cb.Ticks = -0.006:0.003:0.006;

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5 = gca;   % for dx=500 m
ax5.XLim = [xmin,xmax];
ax5.XTick = [xmin:xtick:xmax];

ax5.YLim = [ymin,ymax];
ax5.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(2,2.7,'\(\Delta=500\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=500
%---------------------------------------------------------
end   % end if pn == 6

% END MAKING PLOT OF (-resolved turb. flux of w (FULL PLOT))
%----------------------------------------------------------------------------
% (VII) MAKE PLOT OF Th+Tv-RESOLVED TURB. FLUX OF W (FULL PLOT)
if pn == 7

% FOR dx=125: GET FULL PLOT (NOT DIFFERENCE PLOT)

% To prevent white patches formed on contour plot because of entries in
% thtvturbfxm1r exceeding limits set by clv,
% find those entries and set them as the limit of clv

% initialize new thtvturbfxm1r
thtvturbfxm1rnew = thtvturbfxm1r;

if find(thtvturbfxm1r > clv(end)) ~= 0   % for upper limit of dclv
    thtvturbfxm1rnew(find(thtvturbfxm1r > clv(end))) = clv(end);
end

if find(thtvturbfxm1r < clv(1)) ~= 0   % for lower limit of dclv
    thtvturbfxm1rnew(find(thtvturbfxm1r < clv(1))) = clv(1);
end

subplot(3,1,1)

[C11, h11] = contourf(xh2d1,z2d1,thtvturbfxm1rnew,clv); hold on 
[C21, h21] = contour(xh2d1,z2d1,mywhrm1r,llvp); hold on
[C31, h31] = contour(xh2d1,z2d1,mywhrm1r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h11.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
cb.Ticks = -0.006:0.003:0.006;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\overline{\textrm{T}_h + \textrm{T}_v - \nabla\cdot (\overline{\textrm{\textbf{u}}^{\prime} w^{\prime}})}\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% title(['\(\overline{\textrm{T}_h + \textrm{T}_v - \textrm{resolved turb. flux}}\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% title(['\(\textrm{T}_h + \textrm{T}_v - \textrm{resolved turb. flux} \), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1 = gca;   % for dx=125 m
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(2.5,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(2,2.7,'\(\Delta=125\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% END for dx=125
%---------------------------------------------------------
% FOR dx=250 

subplot(3,1,2)
[C12, h12] = contourf(xh2d2,z2d2,thtvturbfxm2r,clv); hold on 
[C22, h22] = contour(xh2d2,z2d2,mywhrm2r,llvp); hold on
[C32, h32] = contour(xh2d2,z2d2,mywhrm2r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h12.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
cb.Ticks = -0.006:0.003:0.006;
% cb.Ruler.Exponent = -3;
% cb.Ticks = clv(1:2:end);

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2 = gca;   % for dx=250 m
ax2.XLim = [xmin,xmax];
ax2.XTick = [xmin:xtick:xmax];

ax2.YLim = [ymin,ymax];
ax2.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(2,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=250

%---------------------------------------------------------

% FOR dx = 500 m: GET FULL PLOT(NOT DIFFERENCE PLOT)

subplot(3,1,3)

[C15, h15] = contourf(xh2d5,z2d5,thtvturbfxm5r,clv); hold on 
[C25, h25] = contour(xh2d5,z2d5,mywhrm5r,llvp); hold on
[C35, h35] = contour(xh2d5,z2d5,mywhrm5r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h15.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
cb.Ticks = -0.006:0.003:0.006;

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5 = gca;   % for dx=500 m
ax5.XLim = [xmin,xmax];
ax5.XTick = [xmin:xtick:xmax];

ax5.YLim = [ymin,ymax];
ax5.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(2,2.7,'\(\Delta=500\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=500
%---------------------------------------------------------
end   % end if pn == 7

% END MAKING PLOT OF (Th+Tv-resolved turb. flux of w) (FULL PLOT)
%----------------------------------------------------------------------------

% (VIII) MAKE PLOT OF Dw-RESOLVED TURB. FLUX OF W (FULL PLOT)
if pn == 8
    
% To prevent white patches formed on contour plot because of entries in
% dwturbfxm1r exceeding limits set by clv,
% find those entries and set them as the limit of clv

% initialize new dwturbfxm1r
dwturbfxm1rnew = dwturbfxm1r;

if find(dwturbfxm1r > clv(end)) ~= 0   % for upper limit of clv
    dwturbfxm1rnew(find(dwturbfxm1r > clv(end))) = clv(end);
end

if find(dwturbfxm1r < clv(1)) ~= 0   % for lower limit of clv
    dwturbfxm1rnew(find(dwturbfxm1r < clv(1))) = clv(1);
end

% FOR dx=125: GET FULL PLOT

subplot(3,1,1)

[C11, h11] = contourf(xh2d1,z2d1,dwturbfxm1rnew,clv); hold on 
[C21, h21] = contour(xh2d1,z2d1,mywhrm1r,llvp); hold on
[C31, h31] = contour(xh2d1,z2d1,mywhrm1r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h11.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
cb.Ticks = -0.006:0.003:0.006;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\overline{\textrm{D}_w - \nabla\cdot (\overline{\textrm{\textbf{u}}^{\prime} w^{\prime}})}\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% title([ '\(\overline{\textrm{D}_w - \textrm{resolved turb. flux}}\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1 = gca;   % for dx=125 m
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(2.5,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(2,2.7,'\(\Delta=125\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% END for dx=125
%---------------------------------------------------------
% FOR dx=250 

subplot(3,1,2)
[C12, h12] = contourf(xh2d2,z2d2,dwturbfxm2r,clv); hold on 
[C22, h22] = contour(xh2d2,z2d2,mywhrm2r,llvp); hold on
[C32, h32] = contour(xh2d2,z2d2,mywhrm2r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h12.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
% cb.Ruler.Exponent = -3;
cb.Ticks = -0.006:0.003:0.006;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2 = gca;   % for dx=250 m
ax2.XLim = [xmin,xmax];
ax2.XTick = [xmin:xtick:xmax];

ax2.YLim = [ymin,ymax];
ax2.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(2,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=250

%---------------------------------------------------------

% FOR dx = 500 m: GET FULL PLOT

subplot(3,1,3)

[C15, h15] = contourf(xh2d5,z2d5,dwturbfxm5r,clv); hold on 
[C25, h25] = contour(xh2d5,z2d5,mywhrm5r,llvp); hold on
[C35, h35] = contour(xh2d5,z2d5,mywhrm5r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h15.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cb = colorbar;
cb.Ticks = -0.006:0.003:0.006;

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5 = gca;   % for dx=500 m
ax5.XLim = [xmin,xmax];
ax5.XTick = [xmin:xtick:xmax];

ax5.YLim = [ymin,ymax];
ax5.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
text(2.8,3.3,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(2,2.7,'\(\Delta=500\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=500
%---------------------------------------------------------
end   % end if pn == 8

% END MAKING PLOT OF (Dw-resolved turb. flux of w) (FULL PLOT)
%----------------------------------------------------------------------------

% (XI) MAKE PLOT OF PROFILE of 
% (i) Dw-RESOLVED TURB. FLUX OF W 
% (ii) Th+Tv-RESOLVED TURB. FLUX OF W 
%(FULL PLOT, AVG. OVER CENTRAL 1 km on X-Domain)

if pn == 9

% for dx=125
dxh1 = xh1(2:end) - xh1(1:end-1);   % dxh for dx=125 case
xmin1 = find(abs(xh1(xhii1:xhfi1)-(-0.5)) < dxh1(1)/2);   % index of -0.5 on xh1
xmax1 = find(abs(xh1(xhii1:xhfi1)-(0.5)) < dxh1(1)/2);   % index of +0.5 on xh1

% for dx=250
dxh2 = xh2(2:end) - xh2(1:end-1);   % dxh for dx=250 case
xmin2 = find(abs(xh2(xhii2:xhfi2)-(-0.5)) < dxh2(1)/2);   % index of -0.5 on xh2
xmax2 = find(abs(xh2(xhii2:xhfi2)-(0.5)) < dxh2(1)/2);   % index of +0.5 on xh2

% for dx=500
dxh5 = xh5(2:end) - xh5(1:end-1);   % dxh for dx=500 case
xmin5 = max(find(abs(xh5(xhii5:xhfi5)-(-0.5)) < dxh5(1)/1.5));   % index of -0.5 on xh5
xmax5 = min(find(abs(xh5(xhii5:xhfi5)-(0.5)) < dxh5(1)/1.5));   % index of +0.5 on xh5
%-----------------------------------------------------

% get mean y-avg. Dw-RESOLVED TURB. FLUX OF W over central 1 km on x

dwturbfxm1rm = mean(dwturbfxm1r(xmin1:xmax1,:),1);   % for dx=125
dwturbfxm2rm = mean(dwturbfxm2r(xmin2:xmax2,:),1);   % for dx=250
dwturbfxm5rm = mean(dwturbfxm5r(xmin5:xmax5,:),1);   % for dx=250

% get mean y-avg. Th+Tv-RESOLVED TURB. FLUX OF W over central 1 km on x

thtvturbfxm1rm = mean(thtvturbfxm1r(xmin1:xmax1,:),1);   % for dx=125
thtvturbfxm2rm = mean(thtvturbfxm2r(xmin2:xmax2,:),1);   % for dx=250
thtvturbfxm5rm = mean(thtvturbfxm5r(xmin5:xmax5,:),1);   % for dx=250

% add zero level on half-level z to match with size of dwturbfxm2rm matrix
znew = [0; z];

%-----------------------------------------------
% (i) Dw-RESOLVED TURB. FLUX OF W 

subplot(1,2,1)
plot(thtvturbfxm1rm,znew(zii:zfi+1),'--g','LineWidth',1.125,'color',[0 0.5 0]); hold on   % dx=125
plot(thtvturbfxm2rm,znew(zii:zfi+1),'--b','LineWidth',1.125); hold on   % dx=250
plot(thtvturbfxm5rm,znew(zii:zfi+1),'--r','LineWidth',1.125); hold off   % dx=500

grid on;
box on;

lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m');
% lgd = legend('\(\Delta=62.5\) m','\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m','\(\Delta=1000\) m');
set(lgd,'Interpreter','latex','Location','northeast');

xlabel('\big[ms\textsuperscript{-2}\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

title(['\(\overline{\textrm{T}_h + \textrm{T}_v - \nabla\cdot (\overline{\textrm{\textbf{u}}^{\prime} w^{\prime}})}\)'],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = -5e-3;
xmax = 5e-3;
xtick = 5e-3/2;

% get current axes and store it to variable named ax1
ax1 = gca;
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];
ax1.FontSize = 12;

pbr = 2;
% % Make y-axis:x-axis = 3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

%--------------------------------------------------
% (ii) Th+Tv-RESOLVED TURB. FLUX OF W 

subplot(1,2,2)
plot(dwturbfxm1rm,znew(zii:zfi+1),'--g','LineWidth',1.125,'color',[0 0.5 0]); hold on   % dx=125
plot(dwturbfxm2rm,znew(zii:zfi+1),'--b','LineWidth',1.125); hold on   % dx=250
plot(dwturbfxm5rm,znew(zii:zfi+1),'--r','LineWidth',1.125); hold off   % dx=500

grid on;
box on;

lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m');
% lgd = legend('\(\Delta=62.5\) m','\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m','\(\Delta=1000\) m');
set(lgd,'Interpreter','latex','Location','northeast');

xlabel('\big[ms\textsuperscript{-2}\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

title(['\(\overline{\textrm{D}_w - \nabla\cdot (\overline{\textrm{\textbf{u}}^{\prime} w^{\prime}})}\)'],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = -5e-3;
xmax = 5e-3;
xtick = 5e-3/2;

% get current axes and store it to variable named ax1
ax2 = gca;
ax2.XLim = [xmin,xmax];
ax2.XTick = [xmin:xtick:xmax];

ax2.YLim = [ymin,ymax];
ax2.YTick = [ymin:ytick:ymax];
ax2.FontSize = 12;

pbr = 2;
% % Make y-axis:x-axis = 3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

%-----------------------------------------------

% bring the subplots closer together
ax2.Position(1) = ax1.Position(1) + 0.7*ax1.Position(3);

    

end   % end if pn == 9

%-----------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------

% (X) MAKE PLOT OF X-COMP. OF RESOLVED TURB. FLUX OF W
if pn == 10
    
% FOR dx=125: GET DIFFEENCE PLOT

% % interpolate turbfxresxm and myw from dx=125 to dx=250 case
% turbfxresxm1r2 = 0.5.*(turbfxresxm1r(1:2:end-1,:)+ turbfxresxm1r(2:2:end,:));
% mywhrm1r2 = 0.5.*(mywhrm1r(1:2:end-1,:) + mywhrm1r(2:2:end,:));
% 
% % get difference btw. interpolated turbfxresxm1r22 from turbfxresxm2r (mywhrm1r2
% % from mywhrm2r)
% diffturbfxresxm12 = turbfxresxm1r2 - turbfxresxm2r;
% diffmywhrm12 = mywhrm1r2 - mywhrm2r;

%-----------------------
% use interp1 to interpolate grid points from 250-m to 125-m case250
turbfxresxm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));   % y-avg. x-comp. res. turb. flux of w' from dx=250 to dx=125
mywhrm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm1r,2)
turbfxresxm2r1(:,k) = interp1(xh2(xhii2:xhfi2),turbfxresxm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
mywhrm2r1(:,k) = interp1(xh2(xhii2:xhfi2),mywhrm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
end

% get diff. between 125-m and 250-m cases
diffturbfxresxm12 = turbfxresxm1r - turbfxresxm2r1;
diffmywhrm12 = mywhrm1r - mywhrm2r1;   % diff. in w (dx=500-dx=250)

%-----------------------

% subplot(3,1,1)
subplot(3,3,1)

[C11, h11] = contourf(xh2d1(50:110,:),z2d1(50:110,:),diffturbfxresxm12(50:110,:),dclv); hold on 
[C21, h21] = contour(xh2d1(50:110,:),z2d1(50:110,:),diffmywhrm12(50:110,:),dllvp); hold on
[C31, h31] = contour(xh2d1(50:110,:),z2d1(50:110,:),diffmywhrm12(50:110,:),dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h11.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);
% cb = colorbar;

% cb.Ticks = -0.002:0.001:0.002;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs-1);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\overline{-\frac{1}{\rho_{0}}\big[u^{\prime}\frac{\partial}{\partial x}(\rho_{0}w^{\prime}))\big]}\)'],'Interpreter','latex','FontSize',tfs);

% title(['\(\overline{\textrm{T}_h + \textrm{T}_v - \textrm{resolved turb. flux}}\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% title(['\(\textrm{T}_h + \textrm{T}_v - \textrm{resolved turb. flux} \), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1x = gca;   % for dx=125 m
ax1x.XLim = [xmin,xmax];
ax1x.XTick = [xmin:xtick:xmax];

ax1x.YLim = [ymin,ymax];
ax1x.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(2.5,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(2.8,3.3,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % original pos. 

% set panel label at upper left corner and inside the box
text(-2.85,2.6,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % upper left corner

% text(0.85,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(-0.15,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% END for dx=125
%---------------------------------------------------------
% FOR dx=250 

% subplot(3,1,2)
subplot(3,3,4)
[C12, h12] = contourf(xh2d2,z2d2,turbfxresxm2r,clv); hold on 
[C22, h22] = contour(xh2d2,z2d2,mywhrm2r,llvp); hold on
[C32, h32] = contour(xh2d2,z2d2,mywhrm2r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h12.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);

% cb = colorbar;
cb.Ticks = -0.005:0.0025:0.005;
% cb.Ruler.Exponent = -3;
% cb.Ticks = clv(1:2:end);

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs-1);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2x = gca;   % for dx=250 m
ax2x.XLim = [xmin,xmax];
ax2x.XTick = [xmin:xtick:xmax];

ax2x.YLim = [ymin,ymax];
ax2x.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(2.8,3.3,'(b)', 'Interpreter','latex','FontSize',tefs,'FontWeight','bold');   % original position

% set panel label at upper left corner and inside the box
text(-2.85,2.6,'(d)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % upper left corner

% text(2,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(1.6,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=250

%---------------------------------------------------------

% FOR dx = 500 m: GET DIFFERENCE PLOT

% % interpolate turbfxresxm5r and myw from dx=250 to dx=500 case
% turbfxresxm2r5 = 0.5.*(turbfxresxm2r(1:2:end-1,:)+ turbfxresxm2r(2:2:end,:));
% mywhrm2r5 = 0.5.*(mywhrm2r(1:2:end-1,:) + mywhrm2r(2:2:end,:));
% 
% % get difference btw. interpolated turbfxresxm5r from turbfxresxm2r5 (mywhrm5r
% % from mywhrm2r5)
% diffturbfxresxm52 = turbfxresxm5r - turbfxresxm2r5;
% diffmywhrm52 = mywhrm5r - mywhrm2r5;

%-----------------------
% use interp1 to interpolate grid points from 500-m to 250-m case250
turbfxresxm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));   % y-avg. x-comp. res. turb. flux of w' from dx=500 to dx=250
mywhrm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm2r,2)
turbfxresxm5r2(:,k) = interp1(xh5(xhii5:xhfi5),turbfxresxm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
mywhrm5r2(:,k) = interp1(xh5(xhii5:xhfi5),mywhrm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
end

% get diff. between 500-m and 250-m cases
diffturbfxresxm52 = turbfxresxm5r2 - turbfxresxm2r;
diffmywhrm52 = mywhrm5r2 - mywhrm2r;   % diff. in w (dx=500-dx=250)

%-----------------------

% To prevent white patches formed on contour plot because of entries in
% diffturbfxresxm52 exceeding limits set by dclv,
% find those entries and set them as the limit of dclv

% initialize new diffturbfxresxm52
diffturbfxresxm52new = diffturbfxresxm52;

if find(diffturbfxresxm52 > dclv(end)) ~= 0   % for upper limit of dclv
    diffturbfxresxm52new(find(diffturbfxresxm52 > dclv(end))) = dclv(end);
end

if find(diffturbfxresxm52 < dclv(1)) ~= 0   % for lower limit of dclv
    diffturbfxresxm52new(find(diffturbfxresxm52 < dclv(1))) = dclv(1);
end

% subplot(3,1,3)
subplot(3,3,7)

[C15, h15] = contourf(xh2d2,z2d2,diffturbfxresxm52new,dclv); hold on 
[C25, h25] = contour(xh2d2,z2d2,diffmywhrm52,dllvp); hold on
[C35, h35] = contour(xh2d2,z2d2,diffmywhrm52,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h15.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure

% colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);

% cb = colorbar;
% cb.Ticks = -0.002:0.001:0.002;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs-1);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs-1);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5x = gca;   % for dx=500 m
ax5x.XLim = [xmin,xmax];
ax5x.XTick = [xmin:xtick:xmax];

ax5x.YLim = [ymin,ymax];
ax5x.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(2.8,3.3,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % original pos.

% set panel label at upper left corner and inside the box
text(-2.85,2.6,'(g)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % upper left corner

% text(0.85,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(-0.15,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=500
%---------------------------------------------------------
% end   % end if pn == 10

% END MAKING PLOT OF (x-comp. of resolved turb. flux of w)
%----------------------------------------------------------------------------
% (XI) MAKE PLOT OF Y-COMP. OF RESOLVED TURB. FLUX OF W
% if pn == 11

% FOR dx=125: GET DIFFEENCE PLOT

% % interpolate turbfxresxm and myw from dx=125 to dx=250 case
% turbfxresym1r2 = 0.5.*(turbfxresym1r(1:2:end-1,:)+ turbfxresym1r(2:2:end,:));
% mywhrm1r2 = 0.5.*(mywhrm1r(1:2:end-1,:) + mywhrm1r(2:2:end,:));
% 
% % get difference btw. interpolated turbfxresym1r22 from turbfxresym2r (mywhrm1r2
% % from mywhrm2r)
% diffturbfxresym12 = turbfxresym1r2 - turbfxresym2r;
% diffmywhrm12 = mywhrm1r2 - mywhrm2r;

%-----------------------
% use interp1 to interpolate grid points from 250-m to 125-m case250
turbfxresym2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));   % y-avg. y-comp. res. turb. flux of w' from dx=250 to dx=125
mywhrm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm1r,2)
turbfxresym2r1(:,k) = interp1(xh2(xhii2:xhfi2),turbfxresym2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
mywhrm2r1(:,k) = interp1(xh2(xhii2:xhfi2),mywhrm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
end

% get diff. between 125-m and 250-m cases
diffturbfxresym12 = turbfxresym1r - turbfxresym2r1;
diffmywhrm12 = mywhrm1r - mywhrm2r1;   % diff. in w (dx=500-dx=250)

%-----------------------

% subplot(3,1,1)
subplot(3,3,2)

[C11, h11] = contourf(xh2d1(50:110,:),z2d1(50:110,:),diffturbfxresym12(50:110,:),dclv); hold on 
[C21, h21] = contour(xh2d1(50:110,:),z2d1(50:110,:),diffmywhrm12(50:110,:),dllvp); hold on
[C31, h31] = contour(xh2d1(50:110,:),z2d1(50:110,:),diffmywhrm12(50:110,:),dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h11.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);

% cb = colorbar;
% cb.Ticks = -0.002:0.001:0.002;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\overline{-\frac{1}{\rho_{0}}\big[v^{\prime}\frac{\partial}{\partial y}(\rho_{0}w^{\prime}))\big]}\)'],'Interpreter','latex','FontSize',tfs);

% title(['\(\overline{\textrm{T}_h + \textrm{T}_v - \textrm{resolved turb. flux}}\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% title(['\(\textrm{T}_h + \textrm{T}_v - \textrm{resolved turb. flux} \), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1y = gca;   % for dx=125 m
ax1y.XLim = [xmin,xmax];
ax1y.XTick = [xmin:xtick:xmax];

ax1y.YLim = [ymin,ymax];
ax1y.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(2.5,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(2.8,3.3,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % original pos.

% set panel label at upper left corner and inside the box
text(-2.85,2.6,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % upper left corner

% text(0.85,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(-0.15,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% END for dx=125
%---------------------------------------------------------
% FOR dx=250 

% subplot(3,1,2)
subplot(3,3,5)
[C12, h12] = contourf(xh2d2(20:60,:),z2d2(20:60,:),turbfxresym2r(20:60,:),clv); hold on 
[C22, h22] = contour(xh2d2(20:60,:),z2d2(20:60,:),mywhrm2r(20:60,:),llvp); hold on
[C32, h32] = contour(xh2d2(20:60,:),z2d2(20:60,:),mywhrm2r(20:60,:),llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h12.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);

% cb = colorbar;
cb.Ticks = -0.005:0.0025:0.005;
% cb.Ruler.Exponent = -3;
% cb.Ticks = clv(1:2:end);

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2y = gca;   % for dx=250 m
ax2y.XLim = [xmin,xmax];
ax2y.XTick = [xmin:xtick:xmax];

ax2y.YLim = [ymin,ymax];
ax2y.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(2.8,3.3,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % original pos.

% set panel label at upper left corner and inside the box
text(-2.85,2.6,'(e)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % upper left corner

% text(2,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(1.6,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=250

%---------------------------------------------------------

% FOR dx = 500 m: GET DIFFERENCE PLOT

% % interpolate turbfxresym5r and myw from dx=250 to dx=500 case
% turbfxresym2r5 = 0.5.*(turbfxresym2r(1:2:end-1,:)+ turbfxresym2r(2:2:end,:));
% mywhrm2r5 = 0.5.*(mywhrm2r(1:2:end-1,:) + mywhrm2r(2:2:end,:));
% 
% % get difference btw. interpolated turbfxresym5r from turbfxresym2r5 (mywhrm5r
% % from mywhrm2r5)
% diffturbfxresym52 = turbfxresym5r - turbfxresym2r5;
% diffmywhrm52 = mywhrm5r - mywhrm2r5;

%-----------------------
% use interp1 to interpolate grid points from 500-m to 250-m case250
turbfxresym5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));   % y-avg. y-comp. res. turb. flux of w' from dx=500 to dx=250
mywhrm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm2r,2)
turbfxresym5r2(:,k) = interp1(xh5(xhii5:xhfi5),turbfxresym5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
mywhrm5r2(:,k) = interp1(xh5(xhii5:xhfi5),mywhrm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
end

% get diff. between 500-m and 250-m cases
diffturbfxresym52 = turbfxresym5r2 - turbfxresym2r;
diffmywhrm52 = mywhrm5r2 - mywhrm2r;   % diff. in w (dx=500-dx=250)

%-----------------------

% To prevent white patches formed on contour plot because of entries in
% diffturbfxresym52 exceeding limits set by dclv,
% find those entries and set them as the limit of dclv

% initialize new diffturbfxresym52
diffturbfxresym52new = diffturbfxresym52;

if find(diffturbfxresym52 > dclv(end)) ~= 0   % for upper limit of dclv
    diffturbfxresym52new(find(diffturbfxresym52 > dclv(end))) = dclv(end);
end

if find(diffturbfxresym52 < dclv(1)) ~= 0   % for lower limit of dclv
    diffturbfxresym52new(find(diffturbfxresym52 < dclv(1))) = dclv(1);
end

% subplot(3,1,3)
subplot(3,3,8)

[C15, h15] = contourf(xh2d2,z2d2,diffturbfxresym52new,dclv); hold on 
[C25, h25] = contour(xh2d2,z2d2,diffmywhrm52,dllvp); hold on
[C35, h35] = contour(xh2d2,z2d2,diffmywhrm52,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h15.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);

% cb = colorbar;
% cb.Ticks = -0.002:0.001:0.002;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs-1);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5y = gca;   % for dx=500 m
ax5y.XLim = [xmin,xmax];
ax5y.XTick = [xmin:xtick:xmax];

ax5y.YLim = [ymin,ymax];
ax5y.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(2.8,3.3,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % original pos.

% set panel label at upper left corner and inside the box
text(-2.85,2.6,'(h)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % upper left corner

% text(0.85,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(-0.15,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=500
%---------------------------------------------------------
% end   % end if pn == 11

% END MAKING PLOT OF (y-comp. of resolved turb. flux of w)
%----------------------------------------------------------------------------

% (XII) MAKE PLOT OF Z-COMP. OF RESOLVED TURB. FLUX OF W
% if pn == 12

% FOR dx=125: GET DIFFEENCE PLOT

% % interpolate turbfxreszm and myw from dx=125 to dx=250 case
% turbfxreszm1r2 = 0.5.*(turbfxreszm1r(1:2:end-1,:)+ turbfxreszm1r(2:2:end,:));
% mywhrm1r2 = 0.5.*(mywhrm1r(1:2:end-1,:) + mywhrm1r(2:2:end,:));
% 
% % get difference btw. interpolated turbfxreszm1r22 from turbfxreszm2r (mywhrm1r2
% % from mywhrm2r)
% diffturbfxreszm12 = turbfxreszm1r2 - turbfxreszm2r;
% diffmywhrm12 = mywhrm1r2 - mywhrm2r;

%-----------------------
% use interp1 to interpolate grid points from 250-m to 125-m case250
turbfxreszm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));   % y-avg. z-comp. res. turb. flux of w' from dx=250 to dx=125
mywhrm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm1r,2)
turbfxreszm2r1(:,k) = interp1(xh2(xhii2:xhfi2),turbfxreszm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
mywhrm2r1(:,k) = interp1(xh2(xhii2:xhfi2),mywhrm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
end

% get diff. between 125-m and 250-m cases
diffturbfxreszm12 = turbfxreszm1r - turbfxreszm2r1;
diffmywhrm12 = mywhrm1r - mywhrm2r1;   % diff. in w (dx=500-dx=250)

%-----------------------

% subplot(3,1,1)
subplot(3,3,3)

[C11, h11] = contourf(xh2d1(50:110,:),z2d1(50:110,:),diffturbfxreszm12(50:110,:),dclv); hold on 
[C21, h21] = contour(xh2d1(50:110,:),z2d1(50:110,:),diffmywhrm12(50:110,:),dllvp); hold on
[C31, h31] = contour(xh2d1(50:110,:),z2d1(50:110,:),diffmywhrm12(50:110,:),dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h11.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);

cb1z = colorbar;
cb1z.Position(1) = 0.9;

% cb = colorbar;
% cb.Ticks = -0.002:0.001:0.002;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\overline{-\frac{1}{\rho_{0}}\big[w^{\prime}\frac{\partial}{\partial z}(\rho_{0}w^{\prime}))\big]}\)'],'Interpreter','latex','FontSize',tfs);

% title(['\(\overline{\textrm{T}_h + \textrm{T}_v - \textrm{resolved turb. flux}}\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% title(['\(\textrm{T}_h + \textrm{T}_v - \textrm{resolved turb. flux} \), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1z = gca;   % for dx=125 m
ax1z.XLim = [xmin,xmax];
ax1z.XTick = [xmin:xtick:xmax];

ax1z.YLim = [ymin,ymax];
ax1z.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(2.5,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(2.8,3.3,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % original pos.

% set panel label at upper left corner and inside the box
text(-2.85,2.6,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % upper left corner

% text(0.85,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(-0.15,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% END for dx=125
%---------------------------------------------------------
% FOR dx=250 

% subplot(3,1,2)
subplot(3,3,6)

[C12, h12] = contourf(xh2d2(20:60,:),z2d2(20:60,:),turbfxreszm2r(20:60,:),clv); hold on 
[C22, h22] = contour(xh2d2(20:60,:),z2d2(20:60,:),mywhrm2r(20:60,:),llvp); hold on
[C32, h32] = contour(xh2d2(20:60,:),z2d2(20:60,:),mywhrm2r(20:60,:),llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h12.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);

cb2z = colorbar;
cb2z.Position(1) = 0.9;


% cb = colorbar;
cb.Ticks = -0.005:0.0025:0.005;
% cb.Ruler.Exponent = -3;
% cb.Ticks = clv(1:2:end);

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2z = gca;   % for dx=250 m
ax2z.XLim = [xmin,xmax];
ax2z.XTick = [xmin:xtick:xmax];

ax2z.YLim = [ymin,ymax];
ax2z.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(2.8,3.3,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % original pos.

% set panel label at upper left corner and inside the box
text(-2.85,2.6,'(f)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % upper left corner

% text(2,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(1.6,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=250

%---------------------------------------------------------

% FOR dx = 500 m: GET DIFFERENCE PLOT

% % interpolate turbfxreszm5r and myw from dx=250 to dx=500 case
% turbfxreszm2r5 = 0.5.*(turbfxreszm2r(1:2:end-1,:)+ turbfxreszm2r(2:2:end,:));
% mywhrm2r5 = 0.5.*(mywhrm2r(1:2:end-1,:) + mywhrm2r(2:2:end,:));
% 
% % get difference btw. interpolated turbfxreszm5r from turbfxreszm2r5 (mywhrm5r
% % from mywhrm2r5)
% diffturbfxreszm52 = turbfxreszm5r - turbfxreszm2r5;
% diffmywhrm52 = mywhrm5r - mywhrm2r5;

%-----------------------
% use interp1 to interpolate grid points from 500-m to 250-m case250
turbfxreszm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));   % y-avg. z-comp. res. turb. flux of w' from dx=500 to dx=250
mywhrm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm2r,2)
turbfxreszm5r2(:,k) = interp1(xh5(xhii5:xhfi5),turbfxreszm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
mywhrm5r2(:,k) = interp1(xh5(xhii5:xhfi5),mywhrm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
end

% get diff. between 500-m and 250-m cases
diffturbfxreszm52 = turbfxreszm5r2 - turbfxreszm2r;
diffmywhrm52 = mywhrm5r2 - mywhrm2r;   % diff. in w (dx=500-dx=250)

%-----------------------


% To prevent white patches formed on contour plot because of entries in
% diffturbfxreszm52 exceeding limits set by dclv,
% find those entries and set them as the limit of dclv

% initialize new diffturbfxreszm52
diffturbfxreszm52new = diffturbfxreszm52;

if find(diffturbfxreszm52 > dclv(end)) ~= 0   % for upper limit of dclv
    diffturbfxreszm52new(find(diffturbfxreszm52 > dclv(end))) = dclv(end);
end

if find(diffturbfxreszm52 < dclv(1)) ~= 0   % for lower limit of dclv
    diffturbfxreszm52new(find(diffturbfxreszm52 < dclv(1))) = dclv(1);
end

% subplot(3,1,3)
subplot(3,3,9)

[C15, h15] = contourf(xh2d2,z2d2,diffturbfxreszm52new,dclv); hold on 
[C25, h25] = contour(xh2d2,z2d2,diffmywhrm52,dllvp); hold on
[C35, h35] = contour(xh2d2,z2d2,diffmywhrm52,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h15.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);

cb5z = colorbar;
cb5z.Position(1) = 0.9;
% cb = colorbar;
% cb.Ticks = -0.002:0.001:0.002;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs-1);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5z = gca;   % for dx=500 m
ax5z.XLim = [xmin,xmax];
ax5z.XTick = [xmin:xtick:xmax];

ax5z.YLim = [ymin,ymax];
ax5z.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(2.8,3.3,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % original pos.

% set panel label at upper left corner and inside the box
text(-2.85,2.6,'(i)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % upper left corner

% text(0.85,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(-0.15,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=500
%---------------------------------------------------------

% edit subplot positon

% end   % end if pn == 12

% END MAKING PLOT OF (z-comp. of resolved turb. flux of w)
%----------------------------------------------------------------------------

% SET subplots positions

% reduce vertical spacing between them
ax2x.Position(2) = ax5x.Position(2) + ax5x.Position(4) - 0.03;   % for x-comp.
ax1x.Position(2) = ax2x.Position(2) + ax2x.Position(4) - 0.03;

ax2y.Position(2) = ax5y.Position(2) + ax5y.Position(4) - 0.03;   % for y-comp.
ax1y.Position(2) = ax2y.Position(2) + ax2y.Position(4) - 0.03;

ax2z.Position(2) = ax5z.Position(2) + ax5z.Position(4) - 0.03;   % for z-comp.
ax1z.Position(2) = ax2z.Position(2) + ax2z.Position(4) - 0.03;

% reduce horizontal spacing between them 
ax5y.Position(1) = ax5x.Position(1) + ax5x.Position(3) + 0.02;   % for y-comp.
ax2y.Position(1) = ax2x.Position(1) + ax2x.Position(3) + 0.02;   
ax1y.Position(1) = ax1x.Position(1) + ax1x.Position(3) + 0.02; 

ax5z.Position(1) = ax5y.Position(1) + ax5y.Position(3) + 0.02;   % for z-comp.
ax2z.Position(1) = ax2y.Position(1) + ax2y.Position(3) + 0.02;   
ax1z.Position(1) = ax1y.Position(1) + ax1y.Position(3) + 0.02;

% set colorbar width
cb1z.Position(3) = 0.4*cb1z.Position(3);
cb2z.Position(3) = 0.4*cb2z.Position(3);
cb5z.Position(3) = 0.4*cb5z.Position(3);

% reduce spacing between subplots and colorbar
cb1z.Position(1) = ax1z.Position(1) + ax1z.Position(3) + 0.008;
cb2z.Position(1) = ax2z.Position(1) + ax2z.Position(3) + 0.008;
cb5z.Position(1) = ax5z.Position(1) + ax5z.Position(3) + 0.008;

% set colorbar bottom position (level as subplot)
cb1z.Position(2) = ax1z.Position(2) + 0.035;
cb2z.Position(2) = ax2z.Position(2) + 0.035;
cb5z.Position(2) = ax5z.Position(2) + 0.035;

% set height of colorbar (reduce height than original)
cb1z.Position(4) = 0.675*cb1z.Position(4);
cb2z.Position(4) = 0.675*cb2z.Position(4);
cb5z.Position(4) = 0.675*cb5z.Position(4);

end   % end if pn == 10 (for subplots of all comp. of resolved turb. flux of w')

% END MAKING PLOTS
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------

% (XX) MAKE PLOT OF HORI. AND VERT. COMP. OF SUBGRID+RESOLVED TURB. FLUX OF W
if pn == 20

% COMPARING HORI. TURB. FLUX OF W' (Th+resolved)

% FOR dx=125: GET DIFFEENCE PLOT

%-----------------------
if cn == 1
% use interp1 to interpolate grid points from 250-m to 125-m case
hturbmrr2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));   % y-avg. hori. subgrid turb. flux of w' from dx=250 to dx=125
turbfxreshm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));   % y-avg. hori.-comp. res. turb. flux of w' from dx=250 to dx=125
mywhrm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm1r,2)
hturbmrr2r1(:,k) = interp1(xh2(xhii2:xhfi2),hturbmrr2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
turbfxreshm2r1(:,k) = interp1(xh2(xhii2:xhfi2),turbfxreshm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
mywhrm2r1(:,k) = interp1(xh2(xhii2:xhfi2),mywhrm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
end

% get diff. between 125-m and 250-m cases
diffturbfxthresh12 = hturbmrr1r + turbfxreshm1r - (hturbmrr2r1 + turbfxreshm2r1);
diffmywhrm12 = mywhrm1r - mywhrm2r1;   % diff. in w (dx=125-dx=250)

end   % end if cn == 1
%-----------------------
%-----------------------
if cn ~= 1
% use interp1 to interpolate grid points from 250-m to 125-m case
% hturbmrr2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));   % y-avg. hori. subgrid turb. flux of w' from dx=250 to dx=125
turbfxreshm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));   % y-avg. hori.-comp. res. turb. flux of w' from dx=250 to dx=125
mywhrm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm1r,2)
% hturbmrr2r1(:,k) = interp1(xh2(xhii2:xhfi2),hturbmrr2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
turbfxreshm2r1(:,k) = interp1(xh2(xhii2:xhfi2),turbfxreshm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
mywhrm2r1(:,k) = interp1(xh2(xhii2:xhfi2),mywhrm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
end

% get diff. between 125-m and 250-m cases
% diffturbfxthresh12 = hturbmrr1r + turbfxreshm1r - (hturbmrr2r1 + turbfxreshm2r1);
diffturbfxthresh12 = turbfxreshm1r - (turbfxreshm2r1);
diffmywhrm12 = mywhrm1r - mywhrm2r1;   % diff. in w (dx=125-dx=250)
end   % end if cn ~= 1
%-----------------------


% get difference btw. interpolated hori. turb. at 125-m and 250-m (mywhrm1r2
% from mywhrm2r)
% diffturbfxthresh12 = hturbmrr1r2 + turbfxresxm1r2 + turbfxresym1r2 - (hturbmrr2r + turbfxresxm2r + turbfxresym2r);
% diffmywhrm12 = mywhrm1r2 - mywhrm2r;

% subplot(3,1,1)
subplot(3,2,1)

[C11, h11] = contourf(xh2d1,z2d1,diffturbfxthresh12,dclv); hold on 
[C21, h21] = contour(xh2d1,z2d1,diffmywhrm12,dllvp); hold on
[C31, h31] = contour(xh2d1,z2d1,diffmywhrm12,dllvn); hold off

% [C11, h11] = contourf(xh2d2,z2d2,diffturbfxthresh12,dclv); hold on 
% [C21, h21] = contour(xh2d2,z2d2,diffmywhrm12,dllvp); hold on
% [C31, h31] = contour(xh2d2,z2d2,diffmywhrm12,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h11.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);

% cb1h = colorbar;
% cb1h.Position(1) = 0.9;

% cb = colorbar;
% cb.Ticks = -0.002:0.001:0.002;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z} \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(\overline{\textrm{T}_h -\frac{1}{\rho_{0}}\big[w^{\prime}\frac{\partial}{\partial z}(\rho_{0}w^{\prime})\big]}\)'],'Interpreter','latex','FontSize',tfs);
% title(['\(\overline{\textrm{T}_{h,\mathrm{sgs}} - \frac{\partial}{\partial x}(u^{\prime}w^{\prime}) - \frac{\partial}{\partial y}(v^{\prime}w^{\prime})}\)'],'Interpreter','latex','FontSize',tfs);
if cn == 1
    title(['\(\overline{\textrm{T}}_{h,\mathrm{sgs}} - \frac{\partial}{\partial x}(\overline{u^{\prime}w^{\prime}}) - \frac{\partial}{\partial y}(\overline{v^{\prime}w^{\prime}})\)'],'Interpreter','latex','FontSize',tfs);
end
if cn ~= 1
title(['\(- \frac{\partial}{\partial x}(\overline{u^{\prime}w^{\prime}}) - \frac{\partial}{\partial y}(\overline{v^{\prime}w^{\prime}})\)'],'Interpreter','latex','FontSize',tfs);
end

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1h = gca;   % for dx=125 m
ax1h.XLim = [xmin,xmax];
ax1h.XTick = [xmin:xtick:xmax];

ax1h.YLim = [ymin,ymax];
ax1h.YTick = [ymin:ytick:ymax];
ax1h.FontSize = axfs;
ax1h.Title.FontSize = tfs;
ax1h.TickLabelInterpreter = 'latex';

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% set panel label at upper left corner and inside the box
text(-2.90,2.6,'(a)', 'Interpreter', 'latex','FontSize',lbfs,'FontWeight','bold');   % upper left corner

% text(0.85,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(0.7,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% END for dx=125
%---------------------------------------------------------
% FOR dx=250 

% Th+resolved turb. flux of w' at 250-m case

% subplot(3,1,2)
subplot(3,2,3)
if cn == 1
[C12, h12] = contourf(xh2d2,z2d2,hturbmrr2r + turbfxreshm2r,clv); hold on 
end
if cn ~= 1
[C12, h12] = contourf(xh2d2,z2d2,turbfxreshm2r,clv); hold on 
end

[C22, h22] = contour(xh2d2,z2d2,mywhrm2r,llvp); hold on
[C32, h32] = contour(xh2d2,z2d2,mywhrm2r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h12.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);

% cb2h = colorbar;
% cb2h.Position(1) = 0.9;

% cb = colorbar;
% cb.Ticks = -0.005:0.0025:0.005;
% cb.Ruler.Exponent = -3;
% cb.Ticks = clv(1:2:end);

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z} \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2h = gca;   % for dx=250 m
ax2h.XLim = [xmin,xmax];
ax2h.XTick = [xmin:xtick:xmax];

ax2h.YLim = [ymin,ymax];
ax2h.YTick = [ymin:ytick:ymax];

ax2h.FontSize = axfs;
ax2h.TickLabelInterpreter = 'latex';

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% set panel label at upper left corner and inside the box
text(-2.90,2.6,'(c)', 'Interpreter', 'latex','FontSize',lbfs,'FontWeight','bold');   % upper left corner

% text(2,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(1.9,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

if cn == 1 
    text(2.35,0.25,'CTRL', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
end
if cn == 2
    text(2.45,0.25,'DRY', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
end
if cn == 3 
    text(2.27,0.25,'MOIST', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
end

% end for dx=250

%---------------------------------------------------------

% FOR dx = 500 m: GET DIFFERENCE PLOT

%-----------------------
if cn == 1
% use interp1 to interpolate grid points from 500-m to 250-m case
hturbmrr5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));   % y-avg. hori. subgrid turb. flux of w' from dx=500 to dx=250
turbfxreshm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));   % y-avg. hori.-comp. res. turb. flux of w' from dx=500 to dx=250
mywhrm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm2r,2)
hturbmrr5r2(:,k) = interp1(xh5(xhii5:xhfi5),hturbmrr5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
turbfxreshm5r2(:,k) = interp1(xh5(xhii5:xhfi5),turbfxreshm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
mywhrm5r2(:,k) = interp1(xh5(xhii5:xhfi5),mywhrm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
end

% get diff. between 500-m and 250-m cases
diffturbfxthresh52 = hturbmrr5r2 + turbfxreshm5r2 - (hturbmrr2r + turbfxreshm2r);
diffmywhrm52 = mywhrm5r2 - mywhrm2r;   % diff. in w (dx=500-dx=250)
end   % end if cn == 1
%-----------------------
if cn ~= 1
% use interp1 to interpolate grid points from 500-m to 250-m case
% hturbmrr5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));   % y-avg. hori. subgrid turb. flux of w' from dx=500 to dx=250
turbfxreshm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));   % y-avg. hori.-comp. res. turb. flux of w' from dx=500 to dx=250
mywhrm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm2r,2)
% hturbmrr5r2(:,k) = interp1(xh5(xhii5:xhfi5),hturbmrr5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
turbfxreshm5r2(:,k) = interp1(xh5(xhii5:xhfi5),turbfxreshm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
mywhrm5r2(:,k) = interp1(xh5(xhii5:xhfi5),mywhrm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
end

% get diff. between 500-m and 250-m cases
% diffturbfxthresh52 = hturbmrr5r2 + turbfxreshm5r2 - (hturbmrr2r + turbfxreshm2r);
diffturbfxthresh52 = turbfxreshm5r2 - (turbfxreshm2r);
diffmywhrm52 = mywhrm5r2 - mywhrm2r;   % diff. in w (dx=500-dx=250)
end   % end if cn ~= 1
%-----------------------

% get difference btw. interpolated hori. turb. at 500-m and 250-m (mywhrm1r2
% from mywhrm2r)
% diffturbfxthresh52 = hturbmrr5r + turbfxresxm5r + turbfxresym5r - (hturbmrr2r5 + turbfxresxm2r5 + turbfxresym2r5);
% diffmywhrm52 = mywhrm5r - mywhrm2r5;

% To prevent white patches formed on contour plot because of entries in
% diffturbfxreszm52 exceeding limits set by dclv,
% find those entries and set them as the limit of dclv

% initialize new diffturbfxthresh52 
diffturbfxthresh52new = diffturbfxthresh52 ;

if find(diffturbfxthresh52  > dclv(end)) ~= 0   % for upper limit of dclv
    diffturbfxthresh52new(find(diffturbfxthresh52  > dclv(end))) = dclv(end);
end

if find(diffturbfxthresh52  < dclv(1)) ~= 0   % for lower limit of dclv
    diffturbfxthresh52new(find(diffturbfxthresh52  < dclv(1))) = dclv(1);
end

% subplot(3,1,3)
subplot(3,2,5)

[C15, h15] = contourf(xh2d2,z2d2,diffturbfxthresh52new,dclv); hold on 
% [C25, h25] = contour(xh2d5,z2d5,diffmywhrm52,dllvp); hold on
% [C35, h35] = contour(xh2d5,z2d5,diffmywhrm52,dllvn); hold off
[C25, h25] = contour(xh2d2,z2d2,diffmywhrm52,dllvp); hold on
[C35, h35] = contour(xh2d2,z2d2,diffmywhrm52,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h15.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);

% cb5h = colorbar;
% cb5h.Position(1) = 0.9;
% cb = colorbar;
% cb.Ticks = -0.002:0.001:0.002;

% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

ylabel('\textit{z} \big[km\big]','Interpreter','latex','FontSize',lbfs);
xlabel('\textit{x} \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5h = gca;   % for dx=500 m
ax5h.XLim = [xmin,xmax];
ax5h.XTick = [xmin:xtick:xmax];

ax5h.YLim = [ymin,ymax];
ax5h.YTick = [ymin:ytick:ymax];

ax5h.FontSize = axfs;
ax5h.TickLabelInterpreter = 'latex';

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% set panel label at upper left corner and inside the box
text(-2.90,2.6,'(e)', 'Interpreter', 'latex','FontSize',lbfs,'FontWeight','bold');   % upper left corner

% text(0.85,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(0.7,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=500
%---------------------------------------------------------
% END COMPARING HORIZONTAL TURB. FLUX OF W BTW. THE THREE CASES

%------------------------
%------------------------

% COMPARING VERT. TURB. FLUX OF W' (Tv+resolved)

% FOR dx=125: GET DIFFEENCE PLOT

%-----------------------
if cn == 1
% use interp1 to interpolate grid points from 250-m to 125-m case
vturbmrr2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));   % y-avg. vert. subgrid turb. flux of w' from dx=250 to dx=125
turbfxresvm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));   % y-avg. vert.-comp. res. turb. flux of w' from dx=250 to dx=125
mywhrm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm1r,2)
vturbmrr2r1(:,k) = interp1(xh2(xhii2:xhfi2),vturbmrr2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
turbfxresvm2r1(:,k) = interp1(xh2(xhii2:xhfi2),turbfxresvm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
mywhrm2r1(:,k) = interp1(xh2(xhii2:xhfi2),mywhrm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
end

% get diff. between 125-m and 250-m cases
diffturbfxthresv12 = vturbmrr1r + turbfxresvm1r - (vturbmrr2r1 + turbfxresvm2r1);
diffmywhrm12 = mywhrm1r - mywhrm2r1;   % diff. in w (dx=125-dx=250)
end   % end if cn == 1
%-----------------------
if cn ~= 1
% use interp1 to interpolate grid points from 250-m to 125-m case
% vturbmrr2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));   % y-avg. vert. subgrid turb. flux of w' from dx=250 to dx=125
turbfxresvm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));   % y-avg. vert.-comp. res. turb. flux of w' from dx=250 to dx=125
mywhrm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm1r,2)
% vturbmrr2r1(:,k) = interp1(xh2(xhii2:xhfi2),vturbmrr2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
turbfxresvm2r1(:,k) = interp1(xh2(xhii2:xhfi2),turbfxresvm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
mywhrm2r1(:,k) = interp1(xh2(xhii2:xhfi2),mywhrm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
end

% get diff. between 125-m and 250-m cases
% diffturbfxthresv12 = vturbmrr1r + turbfxresvm1r - (vturbmrr2r1 + turbfxresvm2r1);
diffturbfxthresv12 = turbfxresvm1r - (turbfxresvm2r1);
diffmywhrm12 = mywhrm1r - mywhrm2r1;   % diff. in w (dx=125-dx=250)
end   % end if cn ~= 1
%-----------------------


% subplot(3,1,1)
subplot(3,2,2)

[C11, h11] = contourf(xh2d1,z2d1,diffturbfxthresv12,dclv); hold on 
[C21, h21] = contour(xh2d1,z2d1,diffmywhrm12,dllvp); hold on
[C31, h31] = contour(xh2d1,z2d1,diffmywhrm12,dllvn); hold off


grid on 
box on

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h11.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);

cb1v = colorbar;
cb1v.Position(1) = 0.9;
cb1v.Ticks = -0.002:0.001:0.002;
cb1v.FontSize = 11;
cb1v.TickLabelInterpreter = 'latex';

% cb = colorbar;
% cb.Ticks = -0.002:0.001:0.002;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);


% title(['\(\overline{\textrm{T}_{v,\mathrm{sgs}} - \frac{1}{\rho_{0}}\frac{\partial}{\partial z}(\rho_{0}w^{\prime}w^{\prime})}\)'],'Interpreter','latex','FontSize',tfs);
if cn == 1
title(['\(\overline{\textrm{T}}_{v,\mathrm{sgs}} - \frac{1}{\rho_{0}}\frac{\partial}{\partial z}(\rho_{0}\overline{w^{\prime}w^{\prime}})\)'],'Interpreter','latex','FontSize',tfs);
end
if cn ~= 1
title(['\(- \frac{1}{\rho_{0}}\frac{\partial}{\partial z}(\rho_{0}\overline{w^{\prime}w^{\prime}})\)'],'Interpreter','latex','FontSize',tfs);
end

% title(['\(\overline{\textrm{T}_h + \textrm{T}_v - \textrm{resolved turb. flux}}\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1v = gca;   % for dx=125 m
ax1v.XLim = [xmin,xmax];
ax1v.XTick = [xmin:xtick:xmax];

ax1v.YLim = [ymin,ymax];
ax1v.YTick = [ymin:ytick:ymax];

ax1v.FontSize = axfs;
ax1v.Title.FontSize = tfs;
ax1v.TickLabelInterpreter = 'latex';


pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% set panel label at upper left corner and inside the box
text(-2.90,2.6,'(b)', 'Interpreter', 'latex','FontSize',lbfs,'FontWeight','bold');   % upper left corner

% text(0.85,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(0.7,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% END for dx=125
%---------------------------------------------------------
% FOR dx=250 

% Tv+resolved turb. flux of w' at 250-m case

% subplot(3,1,2)
subplot(3,2,4)

if cn == 1
[C12, h12] = contourf(xh2d2,z2d2,vturbmrr2r + turbfxresvm2r,clv); hold on 
end
if cn ~= 1
[C12, h12] = contourf(xh2d2,z2d2,turbfxresvm2r,clv); hold on 
end
[C22, h22] = contour(xh2d2,z2d2,mywhrm2r,llvp); hold on
[C32, h32] = contour(xh2d2,z2d2,mywhrm2r,llvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h12.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [clv(1) clv(end)];   % range of caxis
caxis([cr(1) cr(end)]);

cb2v = colorbar;
cb2v.Position(1) = 0.9;
cb2v.Ticks = -0.005:0.0025:0.005;
cb2v.FontSize = 11;
cb2v.TickLabelInterpreter = 'latex';

% cb = colorbar;
% cb.Ticks = -0.005:0.0025:0.005;
% cb.Ruler.Exponent = -3;
% cb.Ticks = clv(1:2:end);

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2v = gca;   % for dx=250 m
ax2v.XLim = [xmin,xmax];
ax2v.XTick = [xmin:xtick:xmax];

ax2v.YLim = [ymin,ymax];
ax2v.YTick = [ymin:ytick:ymax];

ax2v.FontSize = axfs;
ax2v.TickLabelInterpreter = 'latex';


pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% set panel label at upper left corner and inside the box
text(-2.90,2.6,'(d)', 'Interpreter', 'latex','FontSize',lbfs,'FontWeight','bold');   % upper left corner

% text(2,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(1.90,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

if cn == 1 
    text(2.35,0.25,'CTRL', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
end
if cn == 2
    text(2.45,0.25,'DRY', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
end
if cn == 3 
    text(2.27,0.25,'MOIST', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
end

% end for dx=250

%---------------------------------------------------------

% FOR dx = 500 m: GET DIFFERENCE PLOT

%-----------------------
if cn == 1
% use interp1 to interpolate grid points from 500-m to 250-m case
vturbmrr5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));   % y-avg. vert. subgrid turb. flux of w' from dx=500 to dx=250
turbfxresvm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));   % y-avg. vert.-comp. res. turb. flux of w' from dx=500 to dx=250
mywhrm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm2r,2)
vturbmrr5r2(:,k) = interp1(xh5(xhii5:xhfi5),vturbmrr5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
turbfxresvm5r2(:,k) = interp1(xh5(xhii5:xhfi5),turbfxresvm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
mywhrm5r2(:,k) = interp1(xh5(xhii5:xhfi5),mywhrm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
end

% get diff. between 500-m and 250-m cases
diffturbfxthresv52 = vturbmrr5r2 + turbfxresvm5r2 - (vturbmrr2r + turbfxresvm2r);
diffmywhrm52 = mywhrm5r2 - mywhrm2r;   % diff. in w (dx=500-dx=250)
end   % end if cn == 1
%-----------------------
if cn ~= 1
% use interp1 to interpolate grid points from 500-m to 250-m case
turbfxresvm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));   % y-avg. vert.-comp. res. turb. flux of w' from dx=500 to dx=250
mywhrm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm2r,2)
turbfxresvm5r2(:,k) = interp1(xh5(xhii5:xhfi5),turbfxresvm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
mywhrm5r2(:,k) = interp1(xh5(xhii5:xhfi5),mywhrm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
end

% get diff. between 500-m and 250-m cases
diffturbfxthresv52 = turbfxresvm5r2 - (turbfxresvm2r);
diffmywhrm52 = mywhrm5r2 - mywhrm2r;   % diff. in w (dx=500-dx=250)
end   % end if cn ~= 1
%-----------------------

% To prevent white patches formed on contour plot because of entries in
% diffturbfxthresv52  exceeding limits set by dclv,
% find those entries and set them as the limit of dclv

% initialize new diffturbfxthresv52 
diffturbfxthresv52new = diffturbfxthresv52 ;

if find(diffturbfxthresv52  > dclv(end)) ~= 0   % for upper limit of dclv
    diffturbfxthresv52new(find(diffturbfxthresv52  > dclv(end))) = dclv(end);
end

if find(diffturbfxthresv52  < dclv(1)) ~= 0   % for lower limit of dclv
    diffturbfxthresv52new(find(diffturbfxthresv52  < dclv(1))) = dclv(1);
end

% subplot(3,1,3)
subplot(3,2,6)

[C15, h15] = contourf(xh2d2,z2d2,diffturbfxthresv52new,dclv); hold on 
[C25, h25] = contour(xh2d2,z2d2,diffmywhrm52,dllvp); hold on
[C35, h35] = contour(xh2d2,z2d2,diffmywhrm52,dllvn); hold off

grid on 
box on

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

% set line color in contours of ddx(y-avg. u) to none
h15.LineColor = 'none';   % set dudx contour lines as colorless

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [dclv(1) dclv(end)];   % range of caxis
caxis([cr(1) cr(end)]);

cb5v = colorbar;
cb5v.Position(1) = 0.9;
cb5v.Ticks = -0.002:0.001:0.002;
cb5v.FontSize = 11;
cb5v.TickLabelInterpreter = 'latex';

% cb = colorbar;
% cb.Ticks = -0.002:0.001:0.002;

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
xlabel('\textit{x} \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{u}/dx\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5v = gca;   % for dx=500 m
ax5v.XLim = [xmin,xmax];
ax5v.XTick = [xmin:xtick:xmax];

ax5v.YLim = [ymin,ymax];
ax5v.YTick = [ymin:ytick:ymax];

ax5v.FontSize = axfs;
ax5v.TickLabelInterpreter = 'latex';

pbr = pbra;
% % Make y-axis:x-axis = 1/3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% set panel label at upper left corner and inside the box
text(-2.90,2.6,'(f)', 'Interpreter', 'latex','FontSize',lbfs,'FontWeight','bold');   % upper left corner

% text(0.85,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');
text(0.7,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=500
%---------------------------------------------------------
% END COMPARING VERTICAL TURB. FLUX OF W BTW. THE THREE CASES

%---------------------------------------
% SET subplots positions

% reduce vertical spacing between subplots
ax2h.Position(2) = ax5h.Position(2) + ax5h.Position(4) + 0.0425;
ax1h.Position(2) = ax2h.Position(2) + ax2h.Position(4) + 0.0425;

ax2v.Position(2) = ax5v.Position(2) + ax5v.Position(4) + 0.0425;
ax1v.Position(2) = ax2v.Position(2) + ax2v.Position(4) + 0.0425;

% reduce horizontal spacing between subplots
ax5v.Position(1) = ax5h.Position(1) + ax5h.Position(3) + 0.01;
ax2v.Position(1) = ax2h.Position(1) + ax2h.Position(3) + 0.01;
ax1v.Position(1) = ax1h.Position(1) + ax1h.Position(3) + 0.01;

% set colorbar width
cb1v.Position(3) = 0.45*cb1v.Position(3);
cb2v.Position(3) = 0.45*cb2v.Position(3);
cb5v.Position(3) = 0.45*cb5v.Position(3);

% reduce spacing between subplots and colorbar
cb1v.Position(1) = ax1v.Position(1) + ax1v.Position(3);
cb2v.Position(1) = ax2v.Position(1) + ax2v.Position(3);
cb5v.Position(1) = ax5v.Position(1) + ax5v.Position(3);

% set colorbar bottom position (level as subplot)
cb1v.Position(2) = ax1v.Position(2);
cb2v.Position(2) = ax2v.Position(2);
cb5v.Position(2) = ax5v.Position(2);








end   % end if pn == 20




% save plots
if cn == 1
casestr = 'ctrl';
end
if cn == 2
casestr = 'dry';
end
if cn == 3
casestr = 'moist';
end

if pn == 1
filenameg = sprintf('turbflux_w_2_thtv_%d%d',ti,tf);   % filled contours of Th+Tv and line contours of y-avg. w
end
if pn == 2
filenameg = sprintf('turbflux_w__dw_%d%d',ti,tf);   % filled contours of Dw and line contours of y-avg. w
end
if pn == 3
filenameg = sprintf('turbflux_w_2_returbfx_%d%d',ti,tf);   % filled contours of resolved turb. flux of w and line contours of y-avg. w
end
if pn == 4
filenameg = sprintf('turbflux_w_2_thtvturbfxresw_%d%d_%s',ti,tf,casestr);   % filled contours of (Th+Tv+resolved turb. flux of w) and line contours of y-avg. w
end
if pn == 5
filenameg = sprintf('turbflux_w_2_dwreturbfx_%d%d',ti,tf);   % filled contours of (Dw-resolved turb. flux of w) and line contours of y-avg. w
end
if pn == 6
filenameg = sprintf('turbflux_w_2_returbfx_full_%d%d',ti,tf);   % filled contours of (-resolved turb. flux of w) and line contours of y-avg. w
end
if pn == 7
filenameg = sprintf('turbflux_w_2_thtvreturbfx_full_%d%d',ti,tf);   % filled contours of (Th=Tv-resolved turb. flux of w) and line contours of y-avg. w
end
if pn == 8
filenameg = sprintf('turbflux_w_2_dwreturbfx_full_%d%d',ti,tf);   % filled contours of (Dw-resolved turb. flux of w) and line contours of y-avg. w
end
if pn == 9
filenameg = sprintf('turbflux_w_2_profile_%d%d',ti,tf);   % filled contours of (Dw-resolved turb. flux of w) and line contours of y-avg. w
end

if pn == 10
filenameg = sprintf('turbflux_w_2_turbfxresw_comp_subplot_%d%d',ti,tf);   % filled contours of (each component of resolved turb. flux of w) and line contours of y-avg. w
end
if pn == 20
filenameg = sprintf('turbflux_w_2_thtvturbfxresw_hor_vert_subplot_%d%d',ti,tf);   % filled contours of (hori. and vert. component of total (subgrid+res.) turb. flux of w) and line contours of y-avg. w
end

% saveas(f,fullfile(fnameg,filenameg),'epsc');

end   % end for seq = 6:6