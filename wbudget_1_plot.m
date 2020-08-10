close all
clear all

tic 

% directory for wbudget.mat files
% fname = '/aos/home/mtang/Documents/matrices';
% fname = 'D:\Users\stang33\Desktop';
fname = 'C:\Users\SiuLung\Downloads';  

%--------------------------------------------
% file names for w budegt
filename1 = 'wbudget_1_prcl_dx125_dec03.mat';
filename2 = 'wbudget_1_prcl_dx250_dec03.mat';
filename5 = 'wbudget_1_prcl_dx500_dec03.mat';
%--------------------------------------------
% file names for vertical velocity
filename1w = 'turbflux_w_2_prcl_dx125_dec03.mat';
filename2w = 'turbflux_w_2_prcl_dx250_dec03.mat';
filename5w = 'turbflux_w_2_prcl_dx500_dec03.mat';
%--------------------------------------------

% directory for output graphs
% fnameg = '/aos/home/mtang/Documents/graphs';
% fnameg = 'D:\Users\stang33\Desktop';
fnameg = 'C:\Users\SiuLung\Downloads';

% input range of z interested
zmax = 4;
zmin = 0;

% input plot number interested
% (11) = hori. adv., (12) = vert. adv.
% (21) = hori. diff., (22) = vert. diff., (23) = hori. + vert. diff.
% (31) = hori. turb., (32) = vert. turb.
% (4) = pgf, (5) = rdamp
% (6) = buoy. 
% (7) = r.h.s. terms, (8) = dw/dt
% (9) = profile of r.h.s. and dw/dt
% (10) = total subgrid turb. (hori.+vert.)

pn = 23;

%-------------------------------------------------------------

% % get terms matrcies for w budget
terms1 = load(fullfile(fname,filename1));  % terms array  for dx=125
terms2 = load(fullfile(fname,filename2));  % terms array  for dx=250
terms5 = load(fullfile(fname,filename5));  % terms array  for dx=500 

% % get terms matrcies for vertical velocity
terms1w = load(fullfile(fname,filename1w));  % terms array  for dx=125
terms2w = load(fullfile(fname,filename2w));  % terms array  for dx=250
terms5w = load(fullfile(fname,filename5w));  % terms array  for dx=500
%-------------------------------------------------------------

for seq = 6:6
%    
ti = 0+(seq-1)*1; tf = 1+(seq-1)*1;
% 
% for P = 1:3   % subplot position
% 
fieldname = sprintf('f%d%d',ti,tf);

% get xh and z matrices for all plots 
zf = terms2.terms.(fieldname).zf;  % in km

xh1 = terms1.terms.(fieldname).xh;   % in km for dx=125
xh2 = terms2.terms.(fieldname).xh;   % in km for dx=250
xh5 = terms5.terms.(fieldname).xh;   % in km for dx=500
% xh6 = terms6.terms.(fieldname).xh;   % in km for dx=62.5
% xh10 = terms10.terms.(fieldname).xh;   % in km for dx=1000

% get index of zmax on z (half-level z)
zfi = find(abs(zf-zmax)< 0.1/2);
zii = 1;

% get 2-d xh(xhii:xhfi) and z matrices for making contour plots
xh2d1 = repmat(xh1,[1,zfi-zii+1]);  % dx=125
xh2d2 = repmat(xh2,[1,zfi-zii+1]);  % dx=250
xh2d5 = repmat(xh5,[1,zfi-zii+1]);  % dx=500
% xh2d6 = repmat(xh6,[1,zfi-zii+1]);  % dx=62.5
% xh2d10 = repmat(xh10,[1,zfi-zii+1]);  % dx=1000

zf2d1 = repmat(zf(zii:zfi)',[length(xh1),1]);  % dx=125
zf2d2 = repmat(zf(zii:zfi)',[length(xh2),1]);  % dx=250
zf2d5 = repmat(zf(zii:zfi)',[length(xh5),1]);  % dx=500
% zf2d6 = repmat(zf(zii:zfi)',[length(xh6),1]);  % dx=62.5
% zf2d10 = repmat(zf(zii:zfi)',[length(xh10),1]);  % dx=1000

% get hadvm matrices from arrays
hadvm1 = terms1.terms.(fieldname).hadvm;   % time avg. of y-avg. hori. adv. for dx=125
hadvm2 = terms2.terms.(fieldname).hadvm;   % time avg. of y-avg. hori. adv. for dx=250
hadvm5 = terms5.terms.(fieldname).hadvm;   % time avg. of y-avg. hori. adv. for dx=500
% hadvm10 = terms10.terms.(fieldname).hadvm;   % time avg. of y-avg. hori. adv. for dx=1000
% hadvm6 = terms6.terms.(fieldname).hadvm;   % time avg. of y-avg. hori. adv. for dx=62.5

% reduce hadvm matrices into chosen region
hadvm1r = hadvm1(:,zii:zfi);   % for dx=125
hadvm2r = hadvm2(:,zii:zfi);   % for dx=250
hadvm5r = hadvm5(:,zii:zfi);   % for dx=500
% hadvm10r = hadvm10(:,zii:zfi);   % for dx=1000
% hadvm6r = hadvm6(:,zii:zfi);   % for dx=62.5

% get vadvm matrices from arrays
vadvm1 = terms1.terms.(fieldname).vadvm;   % time avg. of y-avg. vert. adv. for dx=125
vadvm2 = terms2.terms.(fieldname).vadvm;   % time avg. of y-avg. vert. adv. for dx=250
vadvm5 = terms5.terms.(fieldname).vadvm;   % time avg. of y-avg. vert. adv. for dx=500
% vadvm10 = terms10.terms.(fieldname).vadvm;   % time avg. of y-avg. vert. adv. for dx=1000
% vadvm6 = terms6.terms.(fieldname).vadvm;   % time avg. of y-avg. vert. adv. for dx=62.5

% reduce vadvm matrices into chosen region
vadvm1r = vadvm1(:,zii:zfi);   % for dx=125
vadvm2r = vadvm2(:,zii:zfi);   % for dx=250
vadvm5r = vadvm5(:,zii:zfi);   % for dx=500
% vadvm10r = vadvm10(:,zii:zfi);   % for dx=1000
% vadvm6r = vadvm6(:,zii:zfi);   % for dx=62.5

% get hidiffm matrices from arrays
hidiffm1 = terms1.terms.(fieldname).hidiffm;   % time avg. of y-avg. hori. diff. for dx=125
hidiffm2 = terms2.terms.(fieldname).hidiffm;   % time avg. of y-avg. hori. diff. for dx=250
hidiffm5 = terms5.terms.(fieldname).hidiffm;   % time avg. of y-avg. hori. diff. for dx=500
% hidiffm10 = terms10.terms.(fieldname).hidiffm;   % time avg. of y-avg. hori. diff. for dx=1000
% hidiffm6 = terms6.terms.(fieldname).hidiffm;   % time avg. of y-avg. hori. diff. for dx=62.5

% reduce hidiffm matrices into chosen region
hidiffm1r = hidiffm1(:,zii:zfi);   % for dx=125
hidiffm2r = hidiffm2(:,zii:zfi);   % for dx=250
hidiffm5r = hidiffm5(:,zii:zfi);   % for dx=500
% hidiffm10r = hidiffm10(:,zii:zfi);   % for dx=1000
% hidiffm6r = hidiffm6(:,zii:zfi);   % for dx=62.5

% get vidiffm matrices from arrays

vidiffm1 = terms1.terms.(fieldname).vidiffm;   % time avg. of y-avg. vert. diff. for dx=125
vidiffm2 = terms2.terms.(fieldname).vidiffm;   % time avg. of y-avg. vert. diff. for dx=250
vidiffm5 = terms5.terms.(fieldname).vidiffm;   % time avg. of y-avg. vert. diff. for dx=500
% vidiffm10 = terms10.terms.(fieldname).vidiffm;   % time avg. of y-avg. vert. diff. for dx=1000
% vidiffm6 = terms6.terms.(fieldname).vidiffm;   % time avg. of y-avg. vert. diff. for dx=62.5

% reduce vidiffm matrices into chosen region
vidiffm1r = vidiffm1(:,zii:zfi);   % for dx=125
vidiffm2r = vidiffm2(:,zii:zfi);   % for dx=250
vidiffm5r = vidiffm5(:,zii:zfi);   % for dx=500
% vidiffm10r = vidiffm10(:,zii:zfi);   % for dx=1000
% vidiffm6r = vidiffm6(:,zii:zfi);   % for dx=62.5

% get hturbm matrices from arrays

hturbm1 = terms1.terms.(fieldname).hturbm;   % time avg. of y-avg. hori. turb. for dx=125
hturbm2 = terms2.terms.(fieldname).hturbm;   % time avg. of y-avg. hori. turb. for dx=250
hturbm5 = terms5.terms.(fieldname).hturbm;   % time avg. of y-avg. hori. turb. for dx=500
% hturbm10 = terms10.terms.(fieldname).hturbm;   % time avg. of y-avg. hori. turb. for dx=1000
% hturbm6 = terms6.terms.(fieldname).hturbm;   % time avg. of y-avg. hori. turb. for dx=62.5

% reduce hturbm matrices into chosen region
hturbm1r = hturbm1(:,zii:zfi);   % for dx=125
hturbm2r = hturbm2(:,zii:zfi);   % for dx=250
hturbm5r = hturbm5(:,zii:zfi);   % for dx=500
% hturbm10r = hturbm10(:,zii:zfi);   % for dx=1000
% hturbm6r = hturbm6(:,zii:zfi);   % for dx=62.5

% get vturbm matrices from arrays

vturbm1 = terms1.terms.(fieldname).vturbm;   % time avg. of y-avg. vert. turb. for dx=125
vturbm2 = terms2.terms.(fieldname).vturbm;   % time avg. of y-avg. vert. turb. for dx=250
vturbm5 = terms5.terms.(fieldname).vturbm;   % time avg. of y-avg. vert. turb. for dx=500
% vturbm10 = terms10.terms.(fieldname).vturbm;   % time avg. of y-avg. vert. turb. for dx=1000
% vturbm6 = terms6.terms.(fieldname).vturbm;   % time avg. of y-avg. vert. turb. for dx=62.5

% reduce hturbm matrices into chosen region
vturbm1r = vturbm1(:,zii:zfi);   % for dx=125
vturbm2r = vturbm2(:,zii:zfi);   % for dx=250
vturbm5r = vturbm5(:,zii:zfi);   % for dx=500
% vturbm10r = vturbm10(:,zii:zfi);   % for dx=1000
% vturbm6r = vturbm6(:,zii:zfi);   % for dx=62.5

% get pgradm matrices from arrays

pgradm1 = terms1.terms.(fieldname).pgradm;   % time avg. of y-avg. pgf for dx=125
pgradm2 = terms2.terms.(fieldname).pgradm;   % time avg. of y-avg. pgf for dx=250
pgradm5 = terms5.terms.(fieldname).pgradm;   % time avg. of y-avg. pgf for dx=500
% pgradm10 = terms10.terms.(fieldname).pgradm;   % time avg. of y-avg. pgf for dx=1000
% pgradm6 = terms6.terms.(fieldname).pgradm;   % time avg. of y-avg. pgf for dx=62.5

% reduce pgradm matrices into chosen region
pgradm1r = pgradm1(:,zii:zfi);   % for dx=125
pgradm2r = pgradm2(:,zii:zfi);   % for dx=250
pgradm5r = pgradm5(:,zii:zfi);   % for dx=500
% pgradm10r = pgradm10(:,zii:zfi);   % for dx=1000
% pgradm6r = pgradm6(:,zii:zfi);   % for dx=62.5

% get rdampm matrices from arrays

rdampm1 = terms1.terms.(fieldname).rdampm;   % time avg. of y-avg. rdamp for dx=125
rdampm2 = terms2.terms.(fieldname).rdampm;   % time avg. of y-avg. rdamp for dx=250
rdampm5 = terms5.terms.(fieldname).rdampm;   % time avg. of y-avg. rdamp for dx=500
% rdampm10 = terms10.terms.(fieldname).rdampm;   % time avg. of y-avg. rdamp for dx=1000
% rdampm6 = terms6.terms.(fieldname).rdampm;   % time avg. of y-avg. rdamp for dx=62.5

% reduce rdampm matrices into chosen region
rdampm1r = rdampm1(:,zii:zfi);   % for dx=125
rdampm2r = rdampm2(:,zii:zfi);   % for dx=250
rdampm5r = rdampm5(:,zii:zfi);   % for dx=500
% rdampm10r = rdampm10(:,zii:zfi);   % for dx=1000
% rdampm6r = rdampm6(:,zii:zfi);   % for dx=62.5

% get buoym matrices from arrays

buoym1 = terms1.terms.(fieldname).buoym;   % time avg. of y-avg. buoym for dx=125
buoym2 = terms2.terms.(fieldname).buoym;   % time avg. of y-avg. buoym for dx=250
buoym5 = terms5.terms.(fieldname).buoym;   % time avg. of y-avg. buoym for dx=500
% buoym10 = terms10.terms.(fieldname).buoym;   % time avg. of y-avg. buoym for dx=1000
% buoym6 = terms6.terms.(fieldname).buoym;   % time avg. of y-avg. buoym for dx=62.5

% reduce buoym matrices into chosen region
buoym1r = buoym1(:,zii:zfi);   % for dx=125
buoym2r = buoym2(:,zii:zfi);   % for dx=250
buoym5r = buoym5(:,zii:zfi);   % for dx=500
% buoym10r = buoym10(:,zii:zfi);   % for dx=1000
% buoym6r = buoym6(:,zii:zfi);   % for dx=62.5

% get myrhsm matrices from arrays

myrhsm1 = terms1.terms.(fieldname).myrhsm;   % time avg. of y-avg. r.h.s. terms for dx=125
myrhsm2 = terms2.terms.(fieldname).myrhsm;   % time avg. of y-avg. r.h.s. terms for dx=250
myrhsm5 = terms5.terms.(fieldname).myrhsm;   % time avg. of y-avg. r.h.s. terms for dx=500
% myrhsm10 = terms10.terms.(fieldname).myrhsm;   % time avg. of y-avg. r.h.s. terms for dx=1000
% myrhsm6 = terms6.terms.(fieldname).myrhsm;   % time avg. of y-avg. r.h.s. terms for dx=62.5

% reduce myrhsm matrices into chosen region
myrhsm1r = myrhsm1(:,zii:zfi);   % for dx=125
myrhsm2r = myrhsm2(:,zii:zfi);   % for dx=250
myrhsm5r = myrhsm5(:,zii:zfi);   % for dx=500
% myrhsm10r = myrhsm10(:,zii:zfi);   % for dx=1000
% myrhsm6r = myrhsm6(:,zii:zfi);   % for dx=62.5

% get ddtwm matrices from arrays

ddtwm1 = terms1.terms.(fieldname).ddtwm;   % time avg. of y-avg. dwdt for dx=125
ddtwm2 = terms2.terms.(fieldname).ddtwm;   % time avg. of y-avg. dwdt for dx=250
ddtwm5 = terms5.terms.(fieldname).ddtwm;   % time avg. of y-avg. dwdt for dx=500
% ddtwm10 = terms10.terms.(fieldname).ddtwm;   % time avg. of y-avg. dwdt for dx=1000
% ddtwm6 = terms6.terms.(fieldname).ddtwm;   % time avg. of y-avg. dwdt for dx=62.5

% reduce ddtwm matrices into chosen region
ddtwm1r = ddtwm1(:,zii:zfi);   % for dx=125
ddtwm2r = ddtwm2(:,zii:zfi);   % for dx=250
ddtwm5r = ddtwm5(:,zii:zfi);   % for dx=500
% ddtwm10r = ddtwm10(:,zii:zfi);   % for dx=1000
% ddtwm6r = ddtwm6(:,zii:zfi);   % for dx=62.5

%-------------------------------------------
% get y-avg. w matrices from arrays from turbflux_w_2.m
mywhrm1 = terms1w.terms.(fieldname).mywhrm;   % time avg. of y-avg. w for dx=125
mywhrm2 = terms2w.terms.(fieldname).mywhrm;   % time avg. of y-avg. w for dx=250
mywhrm5 = terms5w.terms.(fieldname).mywhrm;   % time avg. of y-avg. w for dx=500

% reduce mywhrm matrices into chosen region
mywhrm1rn = mywhrm1(:,zii:zfi);   % for dx=125
mywhrm2rn = mywhrm2(:,zii:zfi);   % for dx=250
mywhrm5rn = mywhrm5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of mywhrm1rn (to remove white space as z starts at 0.5)
mywhrm1r = [mywhrm1rn(:,1) mywhrm1rn];   % for dx=125
mywhrm2r = [mywhrm2rn(:,1) mywhrm2rn];   % for dx=250
mywhrm5r = [mywhrm5rn(:,1) mywhrm5rn];   % for dx=500

% get xhii and xhfi from each case
z = terms2w.terms.(fieldname).z;  % in km
xhii1 = terms1w.terms.(fieldname).xhii;   % index for xmin for dx=125
xhii2 = terms2w.terms.(fieldname).xhii;   % index for xmin for dx=250
xhii5 = terms5w.terms.(fieldname).xhii;   % index for xmin for dx=500

xhfi1 = terms1w.terms.(fieldname).xhfi;   % index for xmax for dx=125
xhfi2 = terms2w.terms.(fieldname).xhfi;   % index for xmax for dx=250
xhfi5 = terms5w.terms.(fieldname).xhfi;   % index for xmax for dx=500

% get index of zmax on z (half-level z)
zfi = find(abs(z-zmax)< 0.1/2);
zii = 1;
%-------------------------------------------

% set axes font size
axfs = 14;   % axes font size
tfs = 18;   % title font size
tefs = 18;   % text font size (plot number)
lbfs = 18;   % x,y-label font size
cbfs = 12;   % colorbar font size

if pn == 23
axfs = 12;   % axes font size
tfs = 14;   % title font size
tefs = 11;   % text font size (case number, e.g., dx=250 m)
lbfs = 14;   % panel label font size
cbfs = 11;   % colorbar font size
end

v250 = [-0.04875:0.025e-1:0.04875];   % for dx=250
vdiff = 0.5.*v250; % for dx=125-dx=250 and dx=500-dx=250

% % set xmin, xmax, and xticks for all plots
xmina = -3;   % xmin for all plots
xmaxa = 3;   % xmax for all plots
xticka = 1.5;   % xticks for all plots

pbra = 2;   % apsect ratio for plots, y:x ratio

%--------------------------------------------------------------------------------

% % (IA) MAKE HORI. ADV. PLOT

if pn == 11

    v = v250;
    
% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle
subplot(1,3,2);
% [C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),hadvm2r(220:260,:),v);
[C, h1] = contourf(xh2d2(220:260,:),zf2d2(220:260,:),hadvm2r(220:260,:),v);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [v(1) v(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['Hori. Adv., ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

% ax.FontSize = axfs;   % 16 originally

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
%---------------------------------------
% % for dx=62.5 - dx=250
% 
% % interpolate advm from dx=62.5 to dx=250 case
% hadvm6r1 = 0.5.*(hadvm6r(1:2:end-1,:)+ hadvm6r(2:2:end,:));
% hadvm6r2 = 0.5.*(hadvm6r1(1:2:end-1,:)+ hadvm6r1(2:2:end,:));
% 
% % get difference btw. interpolated advm6 from advm2
% diffhadvm62 = hadvm6r2 - hadvm2r;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,1);
% % [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffadvam62new(100:140,:),va12);
% [C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),diffhadvm62(100:140,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=62.5\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% 
% % end for dx=62.5
% %---------------------------------------
% %---------------------------------------

% for dx=125 - dx=250

% interpolate advm from dx=125 to dx=250 case
hadvm1r2 = 0.5.*(hadvm1r(1:2:end-1,:)+ hadvm1r(2:2:end,:));

% get difference btw. interpolated hadvm1 from hadvm2
diffhadvm12 = hadvm1r2 - hadvm2r;

% set contour lines 
va12 = vdiff;

subplot(1,3,1);
% [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffadvam12new(100:140,:),va12);
[C, h1] = contourf(xh2d2(220:260,:),zf2d2(220:260,:),diffhadvm12(220:260,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=125\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=125
%---------------------------------------
%---------------------------------------

% for dx=500 - dx=250

% interpolate hadvm from dx=250 to dx=500 case
hadvmr2r5 = 0.5.*(hadvm2r(1:2:end-1,:)+ hadvm2r(2:2:end,:));

% get difference btw. interpolated hadvm2 from hadvm5
diffhadvm52 = hadvm5r - hadvmr2r5;

% set contour lines 
va12 = vdiff;

subplot(1,3,3);
% [C, h1] = contourf(xh2d5(40:80,:),z2d5(40:80,:),diffadvam52new(40:80,:),va12);
[C, h1] = contourf(xh2d5(100:140,:),zf2d5(100:140,:),diffhadvm52(100:140,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=500\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(d)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end   % end if pn == 11

% end for dx=500
%---------------------------------------
%---------------------------------------
% 
% % for dx=1000 - dx=250
% 
% % interpolate hadvm from dx=250 to dx=1000 case
% hadvmr2r5 = 0.5.*(hadvm2r(1:2:end-1,:)+ hadvm2r(2:2:end,:));
% hadvm2r10 = 0.5.*(hadvmr2r5(1:2:end-1,:)+ hadvmr2r5(2:2:end,:));
% 
% % get difference btw. interpolated hadvm2 from hadvm10
% diffhadvm102 = hadvm10r - hadvm2r10;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,5);
% % [C, h1] = contourf(xh2d10(20:40,:),z2d10(20:40,:),diffadvam102new(20:40,:),va12);
% [C, h1] = contourf(xh2d10(20:40,:),zf2d10(20:40,:),diffhadvm102(20:40,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% % ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=1000\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(e)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

end   % end if pn == 11

% end for dx=1000
%---------------------------------------
% END MAKING HORI. ADV. PLOT
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------

% % (IB) MAKE VERT. ADV. PLOT

if pn == 12

    v = v250;
    
% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle
subplot(1,3,2);
% [C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),vadvm2r(220:260,:),v);
[C, h1] = contourf(xh2d2(220:260,:),zf2d2(220:260,:),vadvm2r(220:260,:),v);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [v(1) v(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['Vert. Adv., ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

% ax.FontSize = axfs;   % 16 originally

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
%---------------------------------------
% % for dx=62.5 - dx=250
% 
% % interpolate advm from dx=62.5 to dx=250 case
% vadvm6r1 = 0.5.*(vadvm6r(1:2:end-1,:)+ vadvm6r(2:2:end,:));
% vadvm6r2 = 0.5.*(vadvm6r1(1:2:end-1,:)+ vadvm6r1(2:2:end,:));
% 
% % get difference btw. interpolated advm6 from advm2
% diffvadvm62 = vadvm6r2 - vadvm2r;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,1);
% % [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffadvam62new(100:140,:),va12);
% [C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),diffvadvm62(100:140,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=62.5\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=62.5
%---------------------------------------
%---------------------------------------

% for dx=125 - dx=250

% interpolate advm from dx=125 to dx=250 case
vadvm1r2 = 0.5.*(vadvm1r(1:2:end-1,:)+ vadvm1r(2:2:end,:));

% get difference btw. interpolated hadvm1 from hadvm2
diffvadvm12 = vadvm1r2 - vadvm2r;

% set contour lines 
va12 = vdiff;

subplot(1,3,1);
% [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffadvam12new(100:140,:),va12);
[C, h1] = contourf(xh2d2(220:280,:),zf2d2(220:280,:),diffvadvm12(220:280,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=125\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=125
%---------------------------------------
%---------------------------------------

% for dx=500 - dx=250

% interpolate hadvm from dx=250 to dx=500 case
vadvmr2r5 = 0.5.*(vadvm2r(1:2:end-1,:)+ vadvm2r(2:2:end,:));

% get difference btw. interpolated hadvm2 from hadvm5
diffvadvm52 = vadvm5r - vadvmr2r5;

% set contour lines 
va12 = vdiff;

subplot(1,3,3);
% [C, h1] = contourf(xh2d5(40:80,:),z2d5(40:80,:),diffadvam52new(40:80,:),va12);
[C, h1] = contourf(xh2d5(100:160,:),zf2d5(100:160,:),diffvadvm52(100:160,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=500\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');  
% text(4.7,3.30,'(d)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end   % end if pn == 11

% end for dx=500
%---------------------------------------
%---------------------------------------
% 
% % for dx=1000 - dx=250
% 
% % interpolate hadvm from dx=250 to dx=1000 case
% vadvmr2r5 = 0.5.*(vadvm2r(1:2:end-1,:)+ vadvm2r(2:2:end,:));
% vadvm2r10 = 0.5.*(vadvmr2r5(1:2:end-1,:)+ vadvmr2r5(2:2:end,:));
% 
% % get difference btw. interpolated hadvm2 from hadvm10
% diffvadvm102 = vadvm10r - vadvm2r10;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,5);
% % [C, h1] = contourf(xh2d10(20:40,:),z2d10(20:40,:),diffadvam102new(20:40,:),va12);
% [C, h1] = contourf(xh2d10(20:40,:),zf2d10(20:40,:),diffvadvm102(20:40,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% % ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=1000\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(e)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

end   % end if pn == 12

% end for dx=1000
%---------------------------------------
% END MAKING VERT. ADV. PLOT
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------

% % (IIA) MAKE HORI. DIFF. PLOT

if pn == 21

    v = v250;
    
% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle
subplot(1,3,2);
% [C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),hadvm2r(220:260,:),v);
[C, h1] = contourf(xh2d2(220:280,:),zf2d2(220:280,:),hidiffm2r(220:280,:),v);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [v(1) v(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['Hori. Diff., ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

% ax.FontSize = axfs;   % 16 originally

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
%---------------------------------------
% % for dx=62.5 - dx=250
% 
% % interpolate hidiffm from dx=62.5 to dx=250 case
% hidiffm6r1 = 0.5.*(hidiffm6r(1:2:end-1,:)+ hidiffm6r(2:2:end,:));
% hidiffm6r2 = 0.5.*(hidiffm6r1(1:2:end-1,:)+ hidiffm6r1(2:2:end,:));
% 
% % get difference btw. interpolated hidiffm6 from hidiffm2
% diffhidiffm62 = hidiffm6r2 - hidiffm2r;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,3,1);
% % [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffhidiffm62new(100:140,:),va12);
% [C, h1] = contourf(xh2d2(220:280,:),zf2d2(220:280,:),diffhidiffm62(220:280,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=62.5\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% 
% % end for dx=62.5
% %---------------------------------------
%---------------------------------------

% for dx=125 - dx=250

% interpolate hidiffm from dx=125 to dx=250 case
hidiffm1r2 = 0.5.*(hidiffm1r(1:2:end-1,:)+ hidiffm1r(2:2:end,:));

% get difference btw. interpolated hidiffm1 from hidiffm2
diffhidiffm12 = hidiffm1r2 - hidiffm2r;

% set contour lines 
va12 = vdiff;

subplot(1,3,1);
% [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffadvam12new(100:140,:),va12);
[C, h1] = contourf(xh2d2(220:280,:),zf2d2(220:280,:),diffhidiffm12(220:280,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=125\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=125
%---------------------------------------
%---------------------------------------

% for dx=500 - dx=250

% interpolate hidiffm from dx=250 to dx=500 case
hidiffm2r5 = 0.5.*(hidiffm2r(1:2:end-1,:)+ hidiffm2r(2:2:end,:));

% get difference btw. interpolated hadvm2 from hadvm5
diffhidiffm52 = hidiffm5r - hidiffm2r5;

% set contour lines 
va12 = vdiff;

subplot(1,3,3);
% [C, h1] = contourf(xh2d5(40:80,:),z2d5(40:80,:),diffadvam52new(40:80,:),va12);
[C, h1] = contourf(xh2d5(100:180,:),zf2d5(100:180,:),diffhidiffm52(100:180,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=500\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(d)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end   % end if pn == 11

% end for dx=500
%---------------------------------------
%---------------------------------------
% 
% for dx=1000 - dx=250
% 
% interpolate hadvm from dx=250 to dx=1000 case
% hidiffm2r5 = 0.5.*(hidiffm2r(1:2:end-1,:)+ hidiffm2r(2:2:end,:));
% hidiffm2r10 = 0.5.*(hidiffm2r5(1:2:end-1,:)+ hidiffm2r5(2:2:end,:));
% 
% get difference btw. interpolated hidiffm2 from hidiffm10
% diffhidiffm102 = hidiffm10r - hidiffm2r10;
% 
% set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,5);
% [C, h1] = contourf(xh2d10(20:40,:),z2d10(20:40,:),diffhidiffm102new(20:40,:),va12);
% [C, h1] = contourf(xh2d10(20:40,:),zf2d10(20:40,:),diffhidiffm102(20:40,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=1000\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(e)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

end   % end if pn == 21

% end for dx=1000
%---------------------------------------
% END MAKING HORI. DIFF. PLOT
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------

% % (IIB) MAKE VERT. DIFF. PLOT

if pn == 22

    v = v250;
    
% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle
subplot(1,3,2);
% [C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),Vidiffm2r(220:260,:),v);
[C, h1] = contourf(xh2d2(200:280,:),zf2d2(200:280,:),vidiffm2r(200:280,:),v);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [v(1) v(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['Vert. Diff., ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

% ax.FontSize = axfs;   % 16 originally

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
%---------------------------------------
% % for dx=62.5 - dx=250
% 
% % interpolate vidiffm from dx=62.5 to dx=250 case
% vidiffm6r1 = 0.5.*(vidiffm6r(1:2:end-1,:)+ vidiffm6r(2:2:end,:));
% vidiffm6r2 = 0.5.*(vidiffm6r1(1:2:end-1,:)+ vidiffm6r1(2:2:end,:));
% 
% % get difference btw. interpolated vidiffm6 from vidiffm2
% diffvidiffm62 = vidiffm6r2 - vidiffm2r;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,1);
% % [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffvidiffm62new(100:140,:),va12);
% [C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),diffvidiffm62(100:140,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=62.5\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% 
% % end for dx=62.5
% %---------------------------------------
%---------------------------------------

% for dx=125 - dx=250

% interpolate vidiffm from dx=125 to dx=250 case
vidiffm1r2 = 0.5.*(vidiffm1r(1:2:end-1,:)+ vidiffm1r(2:2:end,:));

% get difference btw. interpolated vidiffm1 from vidiffm2
diffvidiffm12 = vidiffm1r2 - vidiffm2r;

% set contour lines 
va12 = vdiff;

subplot(1,3,1);
% [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffvidiffm12new(100:140,:),va12);
[C, h1] = contourf(xh2d2(200:280,:),zf2d2(200:280,:),diffvidiffm12(200:280,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=125\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=125
%---------------------------------------
%---------------------------------------

% for dx=500 - dx=250

% interpolate vidiffm from dx=250 to dx=500 case
vidiffm2r5 = 0.5.*(vidiffm2r(1:2:end-1,:)+ vidiffm2r(2:2:end,:));

% get difference btw. interpolated hadvm2 from hadvm5
diffvidiffm52 = vidiffm5r - vidiffm2r5;

% set contour lines 
va12 = vdiff;

subplot(1,3,3);
% [C, h1] = contourf(xh2d5(40:80,:),z2d5(40:80,:),diffvidiffm52new(40:80,:),va12);
[C, h1] = contourf(xh2d5(100:140,:),zf2d5(100:140,:),diffvidiffm52(100:140,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=500\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');  
% text(4.7,3.30,'(d)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end   % end if pn == 11

% end for dx=500
%---------------------------------------
%---------------------------------------
% 
% % for dx=1000 - dx=250
% 
% % interpolate vidiffm from dx=250 to dx=1000 case
% vidiffm2r5 = 0.5.*(vidiffm2r(1:2:end-1,:)+ vidiffm2r(2:2:end,:));
% vidiffm2r10 = 0.5.*(vidiffm2r5(1:2:end-1,:)+ vidiffm2r5(2:2:end,:));
% 
% % get difference btw. interpolated vidiffm2 from vidiffm10
% diffvidiffm102 = vidiffm10r - vidiffm2r10;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,5);
% % [C, h1] = contourf(xh2d10(20:40,:),z2d10(20:40,:),diffvidiffm102new(20:40,:),va12);
% [C, h1] = contourf(xh2d10(20:40,:),zf2d10(20:40,:),diffvidiffm102(20:40,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% % ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=1000\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(e)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

end   % end if pn == 22

% end for dx=1000
%---------------------------------------
% END MAKING VERT. DIFF. PLOT
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------
% % (IIC) MAKE HORI. + VERT. DIFF. PLOT

if pn == 23

%     v = v250;
% v = -7.25e-4:0.5e-4:7.25e-4;
v = -7.125e-4:0.25e-4:7.125e-4;
    
% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
llvp = 0.25:0.5:3.75;   % for positive y-avg. w contours
llvn = -1.75:0.5:-0.25;   % for negative y-avg. w contours   % for negative y-avg. w contours

% for dx=250, place it in the middle
subplot(3,1,2);
[C12, h12] = contourf(xh2d2(220:280,:),zf2d2(220:280,:),hidiffm2r(220:280,:) + vidiffm2r(220:280,:),v); hold on
[C22, h22] = contour(xh2d2(201:280,:),zf2d2(201:280,:),mywhrm2r(:,1:end-1),llvp); hold on
[C32, h32] = contour(xh2d2(201:280,:),zf2d2(201:280,:),mywhrm2r(:,1:end-1),llvn); hold off

% % set line color in contours of ddx(y-avg. u) to none
h12.LineColor = 'none';   % set dudx contour lines as colorless

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [v(1) v(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbar2 = colorbar;
cbar2.Position(1) = 0.9;
cbar2.FontSize = cbfs;
cbar2.Ticks = -7e-4:3.5e-4:7e-4;
cbar2.TickLabelInterpreter = 'latex';

% xlabel('\textit{x} \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z} \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['Numerical Diffusion'],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2 = gca;
ax2.XLim = [xmin,xmax];
ax2.XTick = [xmin:xtick:xmax];

ax2.YLim = [ymin,ymax];
ax2.YTick = [ymin:ytick:ymax];

ax2.FontSize = axfs;
ax2.XLabel.FontSize = lbfs;
ax2.YLabel.FontSize = lbfs;
ax2.TickLabelInterpreter = 'latex';

% pbr = pbra;
pbr = 1/3;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(-2.90,2.6,'(b)', 'Interpreter', 'latex','FontSize',lbfs,'FontWeight','bold');   % upper left corner

text(1.90,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
%---------------------------------------
% % for dx=62.5 - dx=250
% 
% % interpolate hidiffm from dx=62.5 to dx=250 case
% hidiffm6r1 = 0.5.*(hidiffm6r(1:2:end-1,:)+ hidiffm6r(2:2:end,:));
% hidiffm6r2 = 0.5.*(hidiffm6r1(1:2:end-1,:)+ hidiffm6r1(2:2:end,:));
% 
% % get difference btw. interpolated hidiffm6 from hidiffm2
% diffhidiffm62 = hidiffm6r2 - hidiffm2r;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,3,1);
% % [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffhidiffm62new(100:140,:),va12);
% [C, h1] = contourf(xh2d2(220:280,:),zf2d2(220:280,:),diffhidiffm62(220:280,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=62.5\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% 
% % end for dx=62.5
% %---------------------------------------
%---------------------------------------

% for dx=125 - dx=250

dllvp = 0.1:0.2:0.5;   % for diff. in positive y-avg. w contours
dllvn = -0.9:0.2:-0.1;   % for diff. in negative y-avg. w contours

% use interp1 to interpolate grid points from 250-m to 125-m case
hidiffm2r1 = zeros(size(hidiffm1r,1),size(hidiffm1r,2));   % y-avg. hori.+vert. diff. from dx=250 to dx=125
vidiffm2r1 = zeros(size(vidiffm1r,1),size(vidiffm1r,2));   % y-avg. hori.+vert. diff. from dx=250 to dx=125
mywhrm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));  % y-avg. w from dx=250 to dx=125

for k = 1:size(hidiffm1r,2)
hidiffm2r1(:,k) = interp1(xh2,hidiffm2r(:,k),xh1,'linear','extrap');
vidiffm2r1(:,k) = interp1(xh2,vidiffm2r(:,k),xh1,'linear','extrap');
mywhrm2r1(:,k) = interp1(xh2(201:280),mywhrm2r(:,k),xh1(401:560),'linear','extrap');
end

% % get difference btw. interpolated hidiffm2r1 from hidiffm1r
diffnudiffm12 = hidiffm1r + vidiffm1r - (hidiffm2r1 + vidiffm2r1) ;
diffmywhrm12 = mywhrm1r - mywhrm2r1;   % diff. in w (dx=125-dx=250)

% set contour lines 
% va12 = vdiff;
va12 = -10.125e-4:0.25e-4:10.125e-4;

subplot(3,1,1);
[C11, h11] = contourf(xh2d1(450:510,:),zf2d1(450:510,:),diffnudiffm12(450:510,:),va12); hold on
[C21, h21] = contour(xh2d1(401:560,:),zf2d1(401:560,:),diffmywhrm12(:,1:end-1),dllvp); hold on
[C31, h31] = contour(xh2d1(401:560,:),zf2d1(401:560,:),diffmywhrm12(:,1:end-1),dllvn); hold off
% [C, h1] = contourf(xh2d2(220:280,:),zf2d2(220:280,:),diffhidiffm12(220:280,:),va12);

% set line color in contours of ddx(y-avg. u) to none
h11.LineColor = 'none';   % set dudx contour lines as colorless

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbar1 = colorbar;
cbar1.Position(1) = 0.9;
cbar1.FontSize = cbfs;
cbar1.Ticks = -1e-3:0.5e-3:1e-3;
cbar1.TickLabelInterpreter = 'latex';

% xlabel('\textit{x} \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z} \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['Numerical Diffusion'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1 = gca;
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];

ax1.FontSize = axfs;
ax1.XLabel.FontSize = lbfs;
ax1.YLabel.FontSize = lbfs;
ax1.TickLabelInterpreter = 'latex';

% pbr = pbra;
pbr = 1/3;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(-2.90,2.6,'(a)', 'Interpreter', 'latex','FontSize',lbfs,'FontWeight','bold');   % upper left corner

text(0.7,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=125
%---------------------------------------
%---------------------------------------

% for dx=500 - dx=250

% use interp1 to interpolate grid points from 500-m to 250-m case
hidiffm5r2 = zeros(size(hidiffm2r,1),size(hidiffm2r,2));   % y-avg. hori.+vert. diff. from dx=500 to dx=250
vidiffm5r2 = zeros(size(vidiffm2r,1),size(vidiffm2r,2));   % y-avg. hori.+vert. diff. from dx=500 to dx=250
mywhrm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(hidiffm2r,2)
hidiffm5r2(:,k) = interp1(xh5,hidiffm5r(:,k),xh2,'linear','extrap');
vidiffm5r2(:,k) = interp1(xh5,vidiffm5r(:,k),xh2,'linear','extrap');
mywhrm5r2(:,k) = interp1(xh5(101:140),mywhrm5r(:,k),xh2(201:280),'linear','extrap');
end

% % get difference btw. interpolated hidiffm5r2 from hidiffm2r
diffnudiffm52 = hidiffm5r2 + vidiffm5r2 - (hidiffm2r + vidiffm2r) ;
diffmywhrm52 = mywhrm5r2 - mywhrm2r;   % diff. in w (dx=500-dx=250)

% set contour lines 
% va12 = vdiff;
va12 = -10.125e-4:0.25e-4:10.125e-4;

subplot(3,1,3);
% [C, h1] = contourf(xh2d5(100:180,:),zf2d5(100:180,:),diffhidiffm52(100:180,:),va12);
[C1, h15] = contourf(xh2d2(220:260,:),zf2d2(220:260,:),diffnudiffm52(220:260,:),va12); hold on
[C25, h25] = contour(xh2d2(201:280,:),zf2d2(201:280,:),diffmywhrm52(:,1:end-1),dllvp); hold on
[C35, h35] = contour(xh2d2(201:280,:),zf2d2(201:280,:),diffmywhrm52(:,1:end-1),dllvn); hold off

% set line color in contours of ddx(y-avg. u) to none
h15.LineColor = 'none';   % set dudx contour lines as colorless

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbar5 = colorbar;
cbar5.Position(1) = 0.9;
cbar5.FontSize = cbfs;
cbar5.Ticks = -1e-3:0.5e-3:1e-3;
cbar5.TickLabelInterpreter = 'latex';

xlabel('\textit{x} \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z} \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(\Delta=500\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5 = gca;
ax5.XLim = [xmin,xmax];
ax5.XTick = [xmin:xtick:xmax];

ax5.YLim = [ymin,ymax];
ax5.YTick = [ymin:ytick:ymax];

ax5.FontSize = axfs;
ax5.XLabel.FontSize = lbfs;
ax5.YLabel.FontSize = lbfs;
ax5.TickLabelInterpreter = 'latex';

% pbr = pbra;
pbr = 1/3;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(-2.90,2.6,'(c)', 'Interpreter', 'latex','FontSize',lbfs,'FontWeight','bold');   % upper left corner

text(0.7,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=500
%---------------------------------------
%---------------------------------------
% 
% for dx=1000 - dx=250
% 
% interpolate hadvm from dx=250 to dx=1000 case
% hidiffm2r5 = 0.5.*(hidiffm2r(1:2:end-1,:)+ hidiffm2r(2:2:end,:));
% hidiffm2r10 = 0.5.*(hidiffm2r5(1:2:end-1,:)+ hidiffm2r5(2:2:end,:));
% 
% get difference btw. interpolated hidiffm2 from hidiffm10
% diffhidiffm102 = hidiffm10r - hidiffm2r10;
% 
% set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,5);
% [C, h1] = contourf(xh2d10(20:40,:),z2d10(20:40,:),diffhidiffm102new(20:40,:),va12);
% [C, h1] = contourf(xh2d10(20:40,:),zf2d10(20:40,:),diffhidiffm102(20:40,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=1000\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(e)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=1000
%---------------------------------------
%---------------------------------------
% SET subplots positions

% reduce vertical spacing between subplots
ax2.Position(2) = ax5.Position(2) + ax5.Position(4) + 0.0425;
ax1.Position(2) = ax2.Position(2) + ax2.Position(4) + 0.0425;

% set colorbar width
cbar1.Position(3) = 0.45*cbar1.Position(3);
cbar2.Position(3) = 0.45*cbar2.Position(3);
cbar5.Position(3) = 0.45*cbar5.Position(3);

% reduce spacing between subplots and colorbar
cbar1.Position(1) = ax1.Position(1) + ax1.Position(3) - 0.2175;
cbar2.Position(1) = ax2.Position(1) + ax2.Position(3) - 0.2175;
cbar5.Position(1) = ax5.Position(1) + ax5.Position(3) - 0.2175;

% set colorbar bottom position (level as subplot)
cbar1.Position(2) = ax1.Position(2);
cbar2.Position(2) = ax2.Position(2);
cbar5.Position(2) = ax5.Position(2);

% END SETTING SUBPLOT POSITION
%---------------------------------------

end   % end if pn == 23

% END MAKING HORI. + VERT. DIFF. PLOT
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------

% % (IIA) MAKE HORI. TURB. PLOT

if pn == 31

    v = v250;
    
% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle
subplot(1,3,2);
% [C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),hturbm2r(220:260,:),v);
[C, h1] = contourf(xh2d2(200:280,:),zf2d2(200:280,:),hturbm2r(200:280,:),v);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [v(1) v(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['Hori. Turb., ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

% ax.FontSize = axfs;   % 16 originally

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
% %---------------------------------------
% % for dx=62.5 - dx=250
% 
% % interpolate hturbm from dx=62.5 to dx=250 case
% hturbm6r1 = 0.5.*(hturbm6r(1:2:end-1,:)+ hturbm6r(2:2:end,:));
% hturbm6r2 = 0.5.*(hturbm6r1(1:2:end-1,:)+ hturbm6r1(2:2:end,:));
% 
% % get difference btw. interpolated hturbm6 from hturbm2
% diffhturbm62 = hturbm6r2 - hturbm2r;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,1);
% % [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffhidiffm62new(100:140,:),va12);
% [C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),diffhturbm62(100:140,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=62.5\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% text(3.4,3.1,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% % text(4.7,3.30,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% 
% % end for dx=62.5
% %---------------------------------------
%---------------------------------------

% for dx=125 - dx=250

% interpolate hturbm from dx=125 to dx=250 case
hturbm1r2 = 0.5.*(hturbm1r(1:2:end-1,:)+ hturbm1r(2:2:end,:));

% get difference btw. interpolated hturbm1 from hturbm2
diffhturbm12 = hturbm1r2 - hturbm2r;

% set contour lines 
va12 = vdiff;

subplot(1,3,1);
% [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffhturbm12new(100:140,:),va12);
[C, h1] = contourf(xh2d2(200:280,:),zf2d2(200:280,:),diffhturbm12(200:280,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=125\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=125
%---------------------------------------
%---------------------------------------

% for dx=500 - dx=250

% interpolate hturbm from dx=250 to dx=500 case
hturbm2r5 = 0.5.*(hturbm2r(1:2:end-1,:)+ hturbm2r(2:2:end,:));

% get difference btw. interpolated hturbm2 from hturbm5
diffhturbm52 = hturbm5r - hturbm2r5;

% set contour lines 
va12 = vdiff;

subplot(1,3,3);
% [C, h1] = contourf(xh2d5(40:80,:),z2d5(40:80,:),diffadvam52new(40:80,:),va12);
[C, h1] = contourf(xh2d5(100:140,:),zf2d5(100:140,:),diffhturbm52(100:140,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=500\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(d)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end   % end if pn == 11

% end for dx=500
%---------------------------------------
% %---------------------------------------
% 
% % for dx=1000 - dx=250
% 
% % interpolate hturbm from dx=250 to dx=1000 case
% hturbm2r5 = 0.5.*(hturbm2r(1:2:end-1,:)+ hturbm2r(2:2:end,:));
% hturbm2r10 = 0.5.*(hturbm2r5(1:2:end-1,:)+ hturbm2r5(2:2:end,:));
% 
% % get difference btw. interpolated hturbm2 from hturbm10
% diffhturbm102 = hturbm10r - hturbm2r10;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,5);
% % [C, h1] = contourf(xh2d10(20:40,:),z2d10(20:40,:),diffhidiffm102new(20:40,:),va12);
% [C, h1] = contourf(xh2d10(20:40,:),zf2d10(20:40,:),diffhturbm102(20:40,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% % ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=1000\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(e)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

end   % end if pn == 31

% end for dx=1000
%---------------------------------------
% END MAKING HORI. TURB. PLOT
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------

% % (IIIB) MAKE VERT. TURB. PLOT

if pn == 32

    v = v250;
    
% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle
subplot(1,3,2);
% [C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),hturbm2r(220:260,:),v);
[C, h1] = contourf(xh2d2(200:280,:),zf2d2(200:280,:),vturbm2r(200:280,:),v);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [v(1) v(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['Vert. Turb., ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

% ax.FontSize = axfs;   % 16 originally

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
% %---------------------------------------
% % for dx=62.5 - dx=250
% 
% % interpolate vturbm from dx=62.5 to dx=250 case
% vturbm6r1 = 0.5.*(vturbm6r(1:2:end-1,:)+ vturbm6r(2:2:end,:));
% vturbm6r2 = 0.5.*(vturbm6r1(1:2:end-1,:)+ vturbm6r1(2:2:end,:));
% 
% % get difference btw. interpolated vturbm6 from vturbm2
% diffvturbm62 = vturbm6r2 - vturbm2r;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,1);
% % [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffhidiffm62new(100:140,:),va12);
% [C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),diffvturbm62(100:140,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=62.5\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% 
% % end for dx=62.5
% %---------------------------------------
%---------------------------------------

% for dx=125 - dx=250

% interpolate vturbm from dx=125 to dx=250 case
vturbm1r2 = 0.5.*(vturbm1r(1:2:end-1,:)+ vturbm1r(2:2:end,:));

% get difference btw. interpolated vturbm1 from vturbm2
diffvturbm12 = vturbm1r2 - vturbm2r;

% set contour lines 
va12 = vdiff;

subplot(1,3,1);
% [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffvturbm12new(100:140,:),va12);
[C, h1] = contourf(xh2d2(200:280,:),zf2d2(200:280,:),diffvturbm12(200:280,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=125\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=125
%---------------------------------------
%---------------------------------------

% for dx=500 - dx=250

% interpolate vturbm from dx=250 to dx=500 case
vturbm2r5 = 0.5.*(vturbm2r(1:2:end-1,:)+ vturbm2r(2:2:end,:));

% get difference btw. interpolated vturbm2 from vturbm5
diffvturbm52 = vturbm5r - vturbm2r5;

% set contour lines 
va12 = vdiff;

subplot(1,3,3);
% [C, h1] = contourf(xh2d5(40:80,:),z2d5(40:80,:),diffvturbm52new(40:80,:),va12);
[C, h1] = contourf(xh2d5(100:160,:),zf2d5(100:160,:),diffvturbm52(100:160,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=500\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(d)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end   % end if pn == 11

% end for dx=500
%---------------------------------------
%---------------------------------------
% 
% % for dx=1000 - dx=250
% 
% % interpolate vturbm from dx=250 to dx=1000 case
% vturbm2r5 = 0.5.*(vturbm2r(1:2:end-1,:)+ vturbm2r(2:2:end,:));
% vturbm2r10 = 0.5.*(vturbm2r5(1:2:end-1,:)+ vturbm2r5(2:2:end,:));
% 
% % get difference btw. interpolated vturbm2 from vturbm10
% diffvturbm102 = vturbm10r - vturbm2r10;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,5);
% % [C, h1] = contourf(xh2d10(20:40,:),z2d10(20:40,:),diffvturbm102new(20:40,:),va12);
% [C, h1] = contourf(xh2d10(20:40,:),zf2d10(20:40,:),diffvturbm102(20:40,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% % ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=1000\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(e)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

end   % end if pn == 32

% end for dx=1000
%---------------------------------------
% END MAKING VERT. TURB. PLOT
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------

% % (IV) MAKE PGF PLOT

if pn == 4

    v = v250;
    
% figure
f = figure('units','normalized','outerposition',[0 0 1 1])
% f = figure('units','normalized')

%---------------------------------------
% for dx=250, place it in the middle

% To prevent white patches formed on contour plot because of entries in
% pgradm2r exceeding limits set by va12,
% find those entries and set them as the limit of va12

% initialize new pgradm2r
pgradm2rnew = pgradm2r;

if find(pgradm2r > v(end)) ~= 0   % for upper limit of v
    pgradm2rnew(find(pgradm2r > v(end))) = v(end);
end

if find(pgradm2r < v(1)) ~= 0   % for lower limit of v
    pgradm2rnew(find(pgradm2r < v(1))) = v(1);
end

subplot(1,3,2);
% [C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),pgradm2r(220:260,:),v);
% [C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),pgradm2r(100:140,:),v);
[C, h1] = contourf(xh2d2(220:260,:),zf2d2(220:260,:),pgradm2rnew(220:260,:),v);

% set(h1,'LineColor','none')

% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap

blue=[0,0,1];
red=[1,0,0];
white=[1,1,1];
ncolors=11;
nc2=0.5*(ncolors-1);
for nc=1:ncolors
  if (nc <= nc2+1)
    bluered(nc,:)=(nc-1)/nc2*white+(nc2-nc+1)/nc2*blue;
  else
    bluered(nc,:)=(nc-nc2-1)/nc2*red+(ncolors-nc)/nc2*white;
  end
end

colormap(bluered);   %insert chosen colormap

cr = [v(1) v(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbar2 = colorbar;
cbar2.Ticks = -0.048:2*0.008:0.048; 
cbar2.FontSize = cbfs;
cbar2.Ruler.Exponent = -2;
cbar2.Position(3) = 0.7*cbar2.Position(3);
cbar2.TickLabelInterpreter = 'latex';

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
xlabel('\textit{x} \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['pgf, ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);
title(['VPPG'],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2 = gca;
ax2.XLim = [xmin,xmax];
ax2.XTick = [xmin:xtick:xmax];

ax2.YLim = [ymin,ymax];
ax2.YTick = [ymin:ytick:ymax];
ax2.FontSize = axfs;
ax2.XLabel.FontSize = lbfs;
ax2.YLabel.FontSize = lbfs;
ax2.TickLabelInterpreter = 'latex';

% ax.FontSize = axfs;   % 16 originally

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(2.7,3.09,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 

% end for dx=250
%---------------------------------------
%---------------------------------------
% % for dx=62.5 - dx=250
% 
% % interpolate pgradm from dx=62.5 to dx=250 case
% pgradm6r1 = 0.5.*(pgradm6r(1:2:end-1,:)+ pgradm6r(2:2:end,:));
% pgradm6r2 = 0.5.*(pgradm6r1(1:2:end-1,:)+ pgradm6r1(2:2:end,:));
% 
% % get difference btw. interpolated pgradm6 from pgradm2
% diffpgradm62 = pgradm6r2 - pgradm2r;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,1);
% % [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffhidiffm62new(100:140,:),va12);
% [C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),diffpgradm62(100:140,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=62.5\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% 
% % end for dx=62.5
% %---------------------------------------
%---------------------------------------

% for dx=125 - dx=250

% use interp1 to interpolate grid points from 250-m to 125-m case
pgradm2r1 = zeros(size(pgradm1r,1),size(pgradm1r,2));   % y-avg. VPPG from dx=250 to dx=125

for k = 1:size(pgradm1r,2)
pgradm2r1(:,k) = interp1(xh2,pgradm2r(:,k),xh1,'linear','extrap');
end

% get difference btw. interpolated pgradm2r1 from pgradm1
diffpgradm12 = pgradm1r - pgradm2r1;

% set contour lines 
va12 = vdiff;

subplot(1,3,1);
% [C, h1] = contourf(xh2d2(220:260,:),zf2d2(220:260,:),diffpgradm12(220:260,:),va12);
[C, h1] = contourf(xh2d1(450:510,:),zf2d1(450:510,:),diffpgradm12(450:510,:),va12);

% set(h1,'LineColor','none')

% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap

colormap(bluered);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbar1 = colorbar;
cbar1.Ticks = 0.5.*[-0.048:2*0.008:0.048]; 
cbar1.FontSize = cbfs;
cbar1.Ruler.Exponent = -2;
cbar1.Position(3) = 0.7*cbar1.Position(3);
cbar1.TickLabelInterpreter = 'latex';

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
xlabel('\textit{x} \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z} \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=125\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1 = gca;
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];
ax1.FontSize = axfs;
ax1.XLabel.FontSize = lbfs;
ax1.YLabel.FontSize = lbfs;
ax1.TickLabelInterpreter = 'latex';

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(2.7,3.09,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=125
%---------------------------------------
%---------------------------------------

% for dx=500 - dx=250

% use interp1 to interpolate grid points from 500-m to 250-m case
pgradm5r2 = zeros(size(pgradm2r,1),size(pgradm2r,2));   % y-avg. VPPG from dx=500 to dx=250

for k = 1:size(pgradm2r,2)
pgradm5r2(:,k) = interp1(xh5,pgradm5r(:,k),xh2,'linear','extrap');
end

% get difference btw. interpolated pgradm5r2 from pgradm2
diffpgradm52 = pgradm5r2 - pgradm2r;

% set contour lines 
va12 = vdiff;

subplot(1,3,3);
% [C, h1] = contourf(xh2d5(40:80,:),z2d5(40:80,:),diffpgradm52new(40:80,:),va12);
[C, h1] = contourf(xh2d2(220:260,:),zf2d2(220:260,:),diffpgradm52(220:260,:),va12);

% set(h1,'LineColor','none')

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbar5 = colorbar;
cbar5.Ticks = 0.5.*[-0.048:2*0.008:0.048]; 
cbar5.FontSize = cbfs;
cbar5.Ruler.Exponent = -2;
cbar5.Position(3) = 0.7*cbar5.Position(3);
cbar5.TickLabelInterpreter = 'latex';

xlabel('\textit{x} \big[km\big]','Interpreter','latex','FontSize',lbfs);
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=500\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5 = gca;
ax5.XLim = [xmin,xmax];
ax5.XTick = [xmin:xtick:xmax];

ax5.YLim = [ymin,ymax];
ax5.YTick = [ymin:ytick:ymax];
ax5.FontSize = axfs;
ax5.XLabel.FontSize = lbfs;
ax5.YLabel.FontSize = lbfs;
ax5.TickLabelInterpreter = 'latex';

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(2.7,3.09,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');


%-----------------------------------------
% make subplots close to each other

ax1.Position(3) = ax5.Position(3);
ax2.Position(3) = ax5.Position(3);
ax1.Position(4) = ax5.Position(4);
ax2.Position(4) = ax5.Position(4);

% get separation between subplot (a) and (b)
diffpos = ax2.Position(1) - ax1.Position(1) - ax1.Position(3);

% set left position of subplot (b) by reducing separation between the two
% subplots
ax2.Position(1) = ax1.Position(1) + ax1.Position(3) + 0.75*diffpos;

% do the same for subplot (c)
ax5.Position(1) = ax2.Position(1) + ax2.Position(3) + 0.75*diffpos;

cbar1.Position(1) = ax1.Position(1) + 1*ax1.Position(3);
cbar2.Position(1) = ax2.Position(1) + 1*ax2.Position(3);
cbar5.Position(1) = ax5.Position(1) + 1*ax5.Position(3);

%-----------------------------------------

% end   % end if pn == 11

% end for dx=500
%---------------------------------------
% %---------------------------------------
% 
% % for dx=1000 - dx=250
% 
% % interpolate pgradm from dx=250 to dx=1000 case
% pgradm2r5 = 0.5.*(pgradm2r(1:2:end-1,:)+ pgradm2r(2:2:end,:));
% pgradm2r10 = 0.5.*(pgradm2r5(1:2:end-1,:)+ pgradm2r5(2:2:end,:));
% 
% % get difference btw. interpolated pgradm2 from pgradm10
% diffpgradm102 = pgradm10r - pgradm2r10;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,5);
% % [C, h1] = contourf(xh2d10(20:40,:),z2d10(20:40,:),diffpgradm102new(20:40,:),va12);
% [C, h1] = contourf(xh2d10(20:40,:),zf2d10(20:40,:),diffpgradm102(20:40,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% % ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=1000\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(e)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

end   % end if pn == 4

% end for dx=1000
%---------------------------------------
% END MAKING PGF PLOT
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------

% % (V) MAKE RDAMP PLOT

if pn == 5

    v = v250;
    
% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle

% % To prevent white patches formed on contour plot because of entries in
% % pgradm2r exceeding limits set by va12,
% % find those entries and set them as the limit of va12
% 
% % initialize new pgradm2r
% pgradm2rnew = pgradm2r;
% 
% if find(pgradm2r > v(end)) ~= 0   % for upper limit of v
%     pgradm2rnew(find(pgradm2r > v(end))) = v(end);
% end
% 
% if find(pgradm2r < v(1)) ~= 0   % for lower limit of v
%     pgradm2rnew(find(pgradm2r < v(1))) = v(1);
% end

subplot(1,5,3);
% [C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),pgradm2r(220:260,:),v);
[C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),rdampm2r(100:140,:),v);
% [C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),pgradm2rnew(100:140,:),v);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [v(1) v(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['rdamp ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

% ax.FontSize = axfs;   % 16 originally

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
%---------------------------------------
% for dx=62.5 - dx=250

% interpolate rdampm from dx=62.5 to dx=250 case
rdampm6r1 = 0.5.*(rdampm6r(1:2:end-1,:)+ rdampm6r(2:2:end,:));
rdampm6r2 = 0.5.*(rdampm6r1(1:2:end-1,:)+ rdampm6r1(2:2:end,:));

% get difference btw. interpolated rdampm6 from rdampm2
diffrdampm62 = rdampm6r2 - rdampm2r;

% set contour lines 
va12 = vdiff;

subplot(1,5,1);
% [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffhidiffm62new(100:140,:),va12);
[C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),diffrdampm62(100:140,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=62.5\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=62.5
%---------------------------------------
%---------------------------------------

% for dx=125 - dx=250

% interpolate rdampm from dx=125 to dx=250 case
rdampm1r2 = 0.5.*(rdampm1r(1:2:end-1,:)+ rdampm1r(2:2:end,:));

% get difference btw. interpolated rdampm1 from rdampm2
diffrdampm12 = rdampm1r2 - rdampm2r;

% set contour lines 
va12 = vdiff;

subplot(1,5,2);
% [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffvturbm12new(100:140,:),va12);
[C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),diffrdampm12(100:140,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=125\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=125
%---------------------------------------
%---------------------------------------

% for dx=500 - dx=250

% interpolate rdampm from dx=250 to dx=500 case
rdampm2r5 = 0.5.*(rdampm2r(1:2:end-1,:)+ rdampm2r(2:2:end,:));

% get difference btw. interpolated rdampm2 from rdampm5
diffrdampm52 = rdampm5r - rdampm2r5;

% set contour lines 
va12 = vdiff;

subplot(1,5,4);
% [C, h1] = contourf(xh2d5(40:80,:),z2d5(40:80,:),diffpgradm52new(40:80,:),va12);
[C, h1] = contourf(xh2d5(40:80,:),zf2d5(40:80,:),diffrdampm52(40:80,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=500\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(d)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end   % end if pn == 11

% end for dx=500
%---------------------------------------
%---------------------------------------

% for dx=1000 - dx=250

% interpolate rdampm from dx=250 to dx=1000 case
rdampm2r5 = 0.5.*(rdampm2r(1:2:end-1,:)+ rdampm2r(2:2:end,:));
rdampm2r10 = 0.5.*(rdampm2r5(1:2:end-1,:)+ rdampm2r5(2:2:end,:));

% get difference btw. interpolated rdampm2 from rdampm10
diffrdampm102 = rdampm10r - rdampm2r10;

% set contour lines 
va12 = vdiff;

subplot(1,5,5);
% [C, h1] = contourf(xh2d10(20:40,:),z2d10(20:40,:),diffrdampm102new(20:40,:),va12);
[C, h1] = contourf(xh2d10(20:40,:),zf2d10(20:40,:),diffrdampm102(20:40,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=1000\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(e)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

end   % end if pn == 5

% end for dx=1000
%---------------------------------------
% END MAKING RDAMP PLOT
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------

% % (VI) MAKE BUOY. PLOT

if pn == 6

    v = v250;
    
% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle

% % To prevent white patches formed on contour plot because of entries in
% % pgradm2r exceeding limits set by va12,
% % find those entries and set them as the limit of va12
% 
% % initialize new pgradm2r
% pgradm2rnew = pgradm2r;
% 
% if find(pgradm2r > v(end)) ~= 0   % for upper limit of v
%     pgradm2rnew(find(pgradm2r > v(end))) = v(end);
% end
% 
% if find(pgradm2r < v(1)) ~= 0   % for lower limit of v
%     pgradm2rnew(find(pgradm2r < v(1))) = v(1);
% end

subplot(1,3,2);
% [C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),pgradm2r(220:260,:),v);
[C, h1] = contourf(xh2d2(220:260,:),zf2d2(220:260,:),buoym2r(220:260,:),v);
% [C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),pgradm2rnew(100:140,:),v);

% set(h1,'LineColor','none')

% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap

blue=[0,0,1];
red=[1,0,0];
white=[1,1,1];
ncolors=11;
nc2=0.5*(ncolors-1);
for nc=1:ncolors
  if (nc <= nc2+1)
    bluered(nc,:)=(nc-1)/nc2*white+(nc2-nc+1)/nc2*blue;
  else
    bluered(nc,:)=(nc-nc2-1)/nc2*red+(ncolors-nc)/nc2*white;
  end
end

colormap(bluered);   %insert chosen colormap

cr = [v(1) v(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbar2 = colorbar;
cbar2.Ticks = -0.048:2*0.008:0.048; 
cbar2.FontSize = cbfs;
cbar2.Ruler.Exponent = -2;
cbar2.Position(3) = 0.7*cbar2.Position(3);
cbar2.TickLabelInterpreter = 'latex';

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
xlabel('\textit{x} \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['Buoy. ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);
title(['BUOY'],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax2 = gca;
ax2.XLim = [xmin,xmax];
ax2.XTick = [xmin:xtick:xmax];

ax2.YLim = [ymin,ymax];
ax2.YTick = [ymin:ytick:ymax];
ax2.FontSize = axfs;
ax2.XLabel.FontSize = lbfs;
ax2.YLabel.FontSize = lbfs;
ax2.TickLabelInterpreter = 'latex';

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(2.7,3.09,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
%---------------------------------------
% % for dx=62.5 - dx=250
% 
% % interpolate buoym from dx=62.5 to dx=250 case
% buoym6r1 = 0.5.*(buoym6r(1:2:end-1,:)+ buoym6r(2:2:end,:));
% buoym6r2 = 0.5.*(buoym6r1(1:2:end-1,:)+ buoym6r1(2:2:end,:));
% 
% % get difference btw. interpolated buoym6 from buoym2
% diffbuoym62 = buoym6r2 - buoym2r;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,1);
% % [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffbuoym62new(100:140,:),va12);
% [C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),diffbuoym62(100:140,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=62.5\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% 
% % end for dx=62.5
% %---------------------------------------
%---------------------------------------

% for dx=125 - dx=250

%-----------------------------------------
% % use interp1 to interpolate grid points from 250-m to 125-m case
% buoym2r1 = zeros(size(buoym1r,1),size(buoym1r,2));   % y-avg. BUOY from dx=250 to dx=125
% for k = 1:size(buoym1r,2)
% buoym2r1(:,k) = interp1(xh2,buoym2r(:,k),xh1,'linear','extrap');
% end
% 
% % get difference btw. interpolated buoym2r1 from buoym1
% diffbuoym12 = buoym1r - buoym2r1;
%-----------------------------------------

% % interpolate buoym from dx=125 to dx=250 case
buoym1r2 = 0.5.*(buoym1r(1:2:end-1,:)+ buoym1r(2:2:end,:));

% get difference btw. interpolated buoym1 from buoym2
diffbuoym12 = buoym1r2 - buoym2r;
%-----------------------------------------

% set contour lines 
va12 = vdiff;

subplot(1,3,1);
% [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffbuoym12new(100:140,:),va12);
% [C, h1] = contourf(xh2d1(450:510,:),zf2d1(450:510,:),diffbuoym12(450:510,:),va12);
[C, h1] = contourf(xh2d2(220:260,:),zf2d2(220:260,:),diffbuoym12(220:260,:),va12);

% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap

colormap(bluered);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbar1 = colorbar;
cbar1.Ticks = 0.5.*[-0.048:2*0.008:0.048]; 
cbar1.FontSize = cbfs;
cbar1.Ruler.Exponent = -2;
cbar1.Position(3) = 0.7*cbar1.Position(3);
cbar1.TickLabelInterpreter = 'latex';

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
xlabel('\textit{x} \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z} \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=125\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax1 = gca;
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];
ax1.FontSize = axfs;
ax1.XLabel.FontSize = lbfs;
ax1.YLabel.FontSize = lbfs;
ax1.TickLabelInterpreter = 'latex';

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(2.7,3.09,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 

% end for dx=125
%---------------------------------------
%---------------------------------------

% for dx=500 - dx=250

%-----------------------------------------
% % % % use interp1 to interpolate grid points from 500-m to 250-m case
% buoym5r2 = zeros(size(buoym2r,1),size(buoym2r,2));   % y-avg. BUOY from dx=500 to dx=250
% for k = 1:size(buoym2r,2)
% buoym5r2(:,k) = interp1(xh5,buoym5r(:,k),xh2,'linear','extrap');
% end
% 
% % get difference btw. interpolated buoym5r2 from buoym2
% diffbuoym52 = buoym5r2 - buoym2r;  
% %-----------------------------------------

% interpolate buoym from dx=250 to dx=500 case
buoym2r5 = 0.5.*(buoym2r(1:2:end-1,:)+ buoym2r(2:2:end,:));

% get difference btw. interpolated buoym2 from buoym5
diffbuoym52 = buoym5r - buoym2r5;
%-----------------------------------------

% set contour lines 
va12 = vdiff;

subplot(1,3,3);
% [C, h1] = contourf(xh2d2(220:260,:),zf2d2(220:260,:),diffbuoym52(220:260,:),va12);
[C, h1] = contourf(xh2d5(100:140,:),zf2d5(100:140,:),diffbuoym52(100:140,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbar5 = colorbar;
cbar5.Ticks = 0.5.*[-0.048:2*0.008:0.048]; 
cbar5.FontSize = cbfs;
cbar5.Ruler.Exponent = -2;
cbar5.Position(3) = 0.7*cbar5.Position(3);
cbar5.TickLabelInterpreter = 'latex';

xlabel('\textit{x} \big[km\big]','Interpreter','latex','FontSize',lbfs);
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=500\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax5 = gca;
ax5.XLim = [xmin,xmax];
ax5.XTick = [xmin:xtick:xmax];

ax5.YLim = [ymin,ymax];
ax5.YTick = [ymin:ytick:ymax];

ax5.FontSize = axfs;
ax5.XLabel.FontSize = lbfs;
ax5.YLabel.FontSize = lbfs;
ax5.TickLabelInterpreter = 'latex';

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(2.7,3.09,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 


% make subplots close to each other

ax1.Position(3) = ax5.Position(3);
ax2.Position(3) = ax5.Position(3);
ax1.Position(4) = ax5.Position(4);
ax2.Position(4) = ax5.Position(4);


% get separation between subplot (a) and (b)
diffpos = ax2.Position(1) - ax1.Position(1) - ax1.Position(3);

% set left position of subplot (b) by reducing separation between the two
% subplots
ax2.Position(1) = ax1.Position(1) + ax1.Position(3) + 0.75*diffpos;

% do the same for subplot (c)
ax5.Position(1) = ax2.Position(1) + ax2.Position(3) + 0.75*diffpos;

cbar1.Position(1) = ax1.Position(1) + 1*ax1.Position(3);
cbar2.Position(1) = ax2.Position(1) + 1*ax2.Position(3);
cbar5.Position(1) = ax5.Position(1) + 1*ax5.Position(3);

%-----------------------------------------

% end   % end if pn == 11

% end for dx=500
%---------------------------------------
%---------------------------------------
% 
% % for dx=1000 - dx=250
% 
% % interpolate buoym from dx=250 to dx=1000 case
% buoym2r5 = 0.5.*(buoym2r(1:2:end-1,:)+ buoym2r(2:2:end,:));
% buoym2r10 = 0.5.*(buoym2r5(1:2:end-1,:)+ buoym2r5(2:2:end,:));
% 
% % get difference btw. interpolated buoym2 from buoym10
% diffbuoym102 = buoym10r - buoym2r10;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,5);
% % [C, h1] = contourf(xh2d10(20:40,:),z2d10(20:40,:),diffrdampm102new(20:40,:),va12);
% [C, h1] = contourf(xh2d10(20:40,:),zf2d10(20:40,:),diffbuoym102(20:40,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% % ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=1000\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(e)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

end   % end if pn == 6

% end for dx=1000
%---------------------------------------
% END MAKING BUOY. PLOT
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------

% % (VII) MAKE R.H.S. TERMS PLOT

if pn == 7

    v = v250;
    
% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle

% % To prevent white patches formed on contour plot because of entries in
% % pgradm2r exceeding limits set by va12,
% % find those entries and set them as the limit of va12
% 
% % initialize new pgradm2r
% pgradm2rnew = pgradm2r;
% 
% if find(pgradm2r > v(end)) ~= 0   % for upper limit of v
%     pgradm2rnew(find(pgradm2r > v(end))) = v(end);
% end
% 
% if find(pgradm2r < v(1)) ~= 0   % for lower limit of v
%     pgradm2rnew(find(pgradm2r < v(1))) = v(1);
% end

subplot(1,5,3);
% [C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),myrhsm2r(220:260,:),v);
[C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),myrhsm2r(100:140,:),v);
% [C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),pgradm2rnew(100:140,:),v);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [v(1) v(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['R.H.S.' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

% ax.FontSize = axfs;   % 16 originally

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
%---------------------------------------
% for dx=62.5 - dx=250

% interpolate myrhsm from dx=62.5 to dx=250 case
myrhsm6r1 = 0.5.*(myrhsm6r(1:2:end-1,:)+ myrhsm6r(2:2:end,:));
myrhsm6r2 = 0.5.*(myrhsm6r1(1:2:end-1,:)+ myrhsm6r1(2:2:end,:));

% get difference btw. interpolated myrhsm6 from myrhsm2
diffmyrhsm62 = myrhsm6r2 - myrhsm2r;

% set contour lines 
va12 = vdiff;

subplot(1,5,1);
% [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffbuoym62new(100:140,:),va12);
[C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),diffmyrhsm62(100:140,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=62.5\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=62.5
%---------------------------------------
%---------------------------------------

% for dx=125 - dx=250

% interpolate myrhsm from dx=125 to dx=250 case
myrhsm1r2 = 0.5.*(myrhsm1r(1:2:end-1,:)+ myrhsm1r(2:2:end,:));

% get difference btw. interpolated myrhsm1 from myrhsm2
diffmyrhsm12 = myrhsm1r2 - myrhsm2r;

% set contour lines 
va12 = vdiff;

subplot(1,5,2);
% [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffbuoym12new(100:140,:),va12);
[C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),diffmyrhsm12(100:140,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=125\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=125
%---------------------------------------
%---------------------------------------

% for dx=500 - dx=250

% interpolate myrhsm from dx=250 to dx=500 case
myrhsm2r5 = 0.5.*(myrhsm2r(1:2:end-1,:)+ myrhsm2r(2:2:end,:));

% get difference btw. interpolated myrhsm2 from myrhsm5
diffmyrhsm52 = myrhsm5r - myrhsm2r5;

% set contour lines 
va12 = vdiff;

subplot(1,5,4);
% [C, h1] = contourf(xh2d5(40:80,:),z2d5(40:80,:),diffbuoym52new(40:80,:),va12);
[C, h1] = contourf(xh2d5(40:80,:),zf2d5(40:80,:),diffmyrhsm52(40:80,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=500\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(d)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end   % end if pn == 11

% end for dx=500
%---------------------------------------
%---------------------------------------

% for dx=1000 - dx=250

% interpolate myrhsm from dx=250 to dx=1000 case
myrhsm2r5 = 0.5.*(myrhsm2r(1:2:end-1,:)+ myrhsm2r(2:2:end,:));
myrhsm2r10 = 0.5.*(myrhsm2r5(1:2:end-1,:)+ myrhsm2r5(2:2:end,:));

% get difference btw. interpolated myrhsm2 from myrhsm10
diffmyrhsm102 = myrhsm10r - myrhsm2r10;

% set contour lines 
va12 = vdiff;

subplot(1,5,5);
% [C, h1] = contourf(xh2d10(20:40,:),z2d10(20:40,:),diffrdampm102new(20:40,:),va12);
[C, h1] = contourf(xh2d10(20:40,:),zf2d10(20:40,:),diffmyrhsm102(20:40,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=1000\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(e)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

end   % end if pn == 7

% end for dx=1000
%---------------------------------------
% END MAKING R.H.S. TERMS PLOT
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------

% % (VIII) MAKE dw/dt TERMS PLOT

if pn == 8

    v = v250;
    
% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle

% % To prevent white patches formed on contour plot because of entries in
% % pgradm2r exceeding limits set by va12,
% % find those entries and set them as the limit of va12
% 
% % initialize new pgradm2r
% pgradm2rnew = pgradm2r;
% 
% if find(pgradm2r > v(end)) ~= 0   % for upper limit of v
%     pgradm2rnew(find(pgradm2r > v(end))) = v(end);
% end
% 
% if find(pgradm2r < v(1)) ~= 0   % for lower limit of v
%     pgradm2rnew(find(pgradm2r < v(1))) = v(1);
% end

subplot(1,5,3);
% [C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),myrhsm2r(220:260,:),v);
[C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),ddtwm2r(100:140,:),v);
% [C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),pgradm2rnew(100:140,:),v);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [v(1) v(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(d\overline{w}/dt\)' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

% ax.FontSize = axfs;   % 16 originally

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
%---------------------------------------
% for dx=62.5 - dx=250

% interpolate ddtwm from dx=62.5 to dx=250 case
ddtwm6r1 = 0.5.*(ddtwm6r(1:2:end-1,:)+ ddtwm6r(2:2:end,:));
ddtwm6r2 = 0.5.*(ddtwm6r1(1:2:end-1,:)+ ddtwm6r1(2:2:end,:));

% get difference btw. interpolated ddtwm6 from ddtwm2
diffddtwm62 = ddtwm6r2 - ddtwm2r;

% set contour lines 
va12 = vdiff;

subplot(1,5,1);
% [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffbuoym62new(100:140,:),va12);
[C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),diffddtwm62(100:140,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=62.5\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=62.5
%---------------------------------------
%---------------------------------------

% for dx=125 - dx=250

% interpolate ddtwm from dx=125 to dx=250 case
ddtwm1r2 = 0.5.*(ddtwm1r(1:2:end-1,:)+ ddtwm1r(2:2:end,:));

% get difference btw. interpolated ddtwm1 from ddtwm2
diffddtwm12 = ddtwm1r2 - ddtwm2r;

% set contour lines 
va12 = vdiff;

subplot(1,5,2);
% [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffbuoym12new(100:140,:),va12);
[C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),diffddtwm12(100:140,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=125\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=125
%---------------------------------------
%---------------------------------------

% for dx=500 - dx=250

% interpolate ddtwm from dx=250 to dx=500 case
ddtwm2r5 = 0.5.*(ddtwm2r(1:2:end-1,:)+ ddtwm2r(2:2:end,:));

% get difference btw. interpolated ddtwm2 from ddtwm5
diffddtwm52 = ddtwm5r - ddtwm2r5;

% set contour lines 
va12 = vdiff;

subplot(1,5,4);
% [C, h1] = contourf(xh2d5(40:80,:),z2d5(40:80,:),diffbuoym52new(40:80,:),va12);
[C, h1] = contourf(xh2d5(40:80,:),zf2d5(40:80,:),diffddtwm52(40:80,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=500\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(d)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end   % end if pn == 11

% end for dx=500
%---------------------------------------
%---------------------------------------

% for dx=1000 - dx=250

% interpolate ddtwm from dx=250 to dx=1000 case
ddtwm2r5 = 0.5.*(ddtwm2r(1:2:end-1,:)+ ddtwm2r(2:2:end,:));
ddtwm2r10 = 0.5.*(ddtwm2r5(1:2:end-1,:)+ ddtwm2r5(2:2:end,:));

% get difference btw. interpolated ddtwm2 from ddtwm10
diffddtwm102 = ddtwm10r - ddtwm2r10;

% set contour lines 
va12 = vdiff;

subplot(1,5,5);
% [C, h1] = contourf(xh2d10(20:40,:),z2d10(20:40,:),diffddtwm102new(20:40,:),va12);
[C, h1] = contourf(xh2d10(20:40,:),zf2d10(20:40,:),diffddtwm102(20:40,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=1000\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(4.7,3.30,'(e)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

end   % end if pn == 8

% end for dx=1000
%---------------------------------------
% END MAKING  dw/dt PLOT
%--------------------------------------------------------------------------------

% (XI) MAKE PROFILE OF R.H.S. and dw/dt plot

if pn == 9
    
% get indices of -0.5 and +0.5 km on xh 
dxh2 = xh2(2:end) - xh2(1:end-1);   % dxh for dx=250 case
xmin2 = find(abs(xh2-(-0.5)) < dxh2(1)/2);   % index of -0.5 on xh2
xmax2 = find(abs(xh2-(0.5)) < dxh2(1)/2);   % index of +0.5 on xh2


% dxh10 = xh10(2:end) - xh10(1:end-1);   % dxh for dx=1000 case
% xmin10 = find(abs(xh10-(-0.5)) < dxh10(1)/2);   % index of -0.5 on xh10
% xmax10 = find(abs(xh10-(0.5)) < dxh10(1)/2);   % index of +0.5 on xh10

% get hori. mean of y-avg. r.h.s. terms over central 1 km in x and along y
myrhsm2rm = mean(myrhsm2r(xmin2:xmax2,:),1);   % for dx=250
% myrhsm10rm = mean(myrhsm10r(xmin10:xmax10,:),1);   % for dx=1000
% myrhsm10rmturb = mean(myrhsm10rturb(xmin10:xmax10,:),1);   % for dx=1000 with doturbdiag

% get total of y-avg. r.h.s. terms over central 1 km in x and along y
myrhsm2rt = sum(myrhsm2r(xmin2:xmax2,:),1);

% get hori. mean of y-avg. dw/dt over central 1 km in x and along y
ddtwm2rm = mean(ddtwm2r(xmin2:xmax2,:),1);   % for dx=250
% ddtwm10rm = mean(ddtwm10r(xmin10:xmax10,:),1);   % for dx=1000
% ddtwm10rmturb = mean(ddtwm10rturb(xmin10:xmax10,:),1);   % for dx=1000 with doturbdiag

% get total of y-avg. dw/dt over central 1 km in x and along y
ddtwm2rt = sum(ddtwm2r(xmin2:xmax2,:),1);

f = figure('units','normalized','outerposition',[0 0 1 1])

% subplot(1,5,3);
% subplot(1,2,1);
% plot(myrhsm2rm,zf(zii:zfi)); hold on
% plot(ddtwm2rm,zf(zii:zfi)); hold off

subplot(1,2,1);
plot(myrhsm10rm,zf(zii:zfi)); hold on
plot(ddtwm10rm,zf(zii:zfi)); hold off


subplot(1,2,2);
plot(myrhsm10rmturb,zf(zii:zfi)); hold on
plot(ddtwm10rmturb,zf(zii:zfi)); hold off


% subplot(1,2,2);
% plot(myrhsm2rt,zf(zii:zfi)); hold on
% plot(ddtwm2rt,zf(zii:zfi)); hold off

grid on;
box on;

lgd = legend('r.h.s.','\(d\overline{w}/dt\)');
set(lgd,'Interpreter','latex','Location','northeast');

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title([num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;

% xmin = -4e-5;
% xmax = 4e-5;
% xtick = 4e-5/2;

% xmin = -5e-5;
% xmax = 5e-5;
% xtick = 5e-5/2;

xmin = -9.5e-4;
xmax = 9.5e-4;
xtick = 9.5e-4/2;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

    
end   % end if pn == 9


%----------------------------------------------------------------------------

% % (X) MAKE SUBGRID TURB. (HORI.+VERT.) PLOT

if pn == 10

    v = 0.1.*v250;
    
% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle
subplot(1,3,2);
% [C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),hturbm2r(220:260,:),v);
[C, h1] = contourf(xh2d2(220:260,:),zf2d2(220:260,:),hturbm2r(220:260,:)+vturbm2r(220:260,:),v);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [v(1) v(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['Subgrid Turb., ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
    ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

% ax.FontSize = axfs;   % 16 originally

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
% %---------------------------------------
% % for dx=62.5 - dx=250
% 
% % interpolate hturbm from dx=62.5 to dx=250 case
% hturbm6r1 = 0.5.*(hturbm6r(1:2:end-1,:)+ hturbm6r(2:2:end,:));
% hturbm6r2 = 0.5.*(hturbm6r1(1:2:end-1,:)+ hturbm6r1(2:2:end,:));
% 
% % get difference btw. interpolated hturbm6 from hturbm2
% diffhturbm62 = hturbm6r2 - hturbm2r;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,1);
% % [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffhidiffm62new(100:140,:),va12);
% [C, h1] = contourf(xh2d2(100:140,:),zf2d2(100:140,:),diffhturbm62(100:140,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=62.5\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% 
% % end for dx=62.5
% %---------------------------------------
%---------------------------------------

% for dx=125 - dx=250

% interpolate hturbm from dx=125 to dx=250 case
hturbm1r2 = 0.5.*(hturbm1r(1:2:end-1,:)+ hturbm1r(2:2:end,:));
vturbm1r2 = 0.5.*(vturbm1r(1:2:end-1,:)+ vturbm1r(2:2:end,:));

% get difference btw. interpolated hturbm1+vturbm1 from hturbm2+vturbm2
diffsgturbm12 = hturbm1r2+vturbm1r2 - hturbm2r - vturbm2r;

% set contour lines 
va12 = 0.1.*vdiff;

subplot(1,3,1);
% [C, h1] = contourf(xh2d2(100:140,:),z2d2(100:140,:),diffhturbm12new(100:140,:),va12);
[C, h1] = contourf(xh2d2(220:260,:),zf2d2(220:260,:),diffsgturbm12(220:260,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=125\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=125
%---------------------------------------
%---------------------------------------

% for dx=500 - dx=250

% interpolate hturbm+vturbm from dx=250 to dx=500 case
hturbm2r5 = 0.5.*(hturbm2r(1:2:end-1,:)+ hturbm2r(2:2:end,:));
vturbm2r5 = 0.5.*(vturbm2r(1:2:end-1,:)+ vturbm2r(2:2:end,:));

% get difference btw. interpolated hturbm2+vturbm2 from hturbm5+vturbm5
diffsgturbm52 = hturbm5r+vturbm5r - hturbm2r5-vturbm2r5;

% set contour lines 
va12 = 0.1.*vdiff;

subplot(1,3,3);
% [C, h1] = contourf(xh2d5(40:80,:),z2d5(40:80,:),diffadvam52new(40:80,:),va12);
[C, h1] = contourf(xh2d5(100:140,:),zf2d5(100:140,:),diffsgturbm52(100:140,:),va12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [va12(1) va12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\Delta=500\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = xmina;
xmax = xmaxa;
xtick = xticka;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.1,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(d)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end   % end if pn == 11

% end for dx=500
%---------------------------------------
%---------------------------------------

% % for dx=1000 - dx=250
% 
% % interpolate hturbm from dx=250 to dx=1000 case
% hturbm2r5 = 0.5.*(hturbm2r(1:2:end-1,:)+ hturbm2r(2:2:end,:));
% hturbm2r10 = 0.5.*(hturbm2r5(1:2:end-1,:)+ hturbm2r5(2:2:end,:));
% 
% % get difference btw. interpolated hturbm2 from hturbm10
% diffhturbm102 = hturbm10r - hturbm2r10;
% 
% % set contour lines 
% va12 = vdiff;
% 
% subplot(1,5,5);
% % [C, h1] = contourf(xh2d10(20:40,:),z2d10(20:40,:),diffhidiffm102new(20:40,:),va12);
% [C, h1] = contourf(xh2d10(20:40,:),zf2d10(20:40,:),diffhturbm102(20:40,:),va12);
% 
% cm = load('cmap_bous.mat');   %load edited color map
% cmap = cm.cmap;          %access variable from structure
% colormap(cmap);   %insert chosen colormap
% 
% cr = [va12(1) va12(end)];   % range of caxis
% caxis([cr(1) cr(end)]);
% colorbar
% 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% % ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);
% 
% timehi = sprintf('%02d',ti);
% timemini = sprintf('%02d',0);
% timehf = sprintf('%02d',tf);
% timeminf = sprintf('%02d',0);
% 
% title(['\(\Delta=1000\) m \(-\) \(\Delta=250\) m'],'Interpreter','latex','FontSize',tfs);
% % format axes
% ymin = 0;
% ymax = 3;  
% ytick = 1;
% 
% xmin = xmina;
% xmax = xmaxa;
% xtick = xticka;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% pbr = pbra;
% % % Make y-axis:x-axis = 2
% pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% % text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
% text(4.7,3.30,'(e)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

end   % end if pn == 10

% end for dx=1000
%---------------------------------------
% END MAKING SUBGRID TURB. PLOT
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------



% save plots
if pn == 11 
filenameg = sprintf('wbudget_1_hadvm_%d%d',ti,tf);   % hori. adv. 
end
if pn == 12 
filenameg = sprintf('wbudget_1_vadvm_%d%d',ti,tf);   % vert. adv.
end
if pn == 21 
filenameg = sprintf('wbudget_1_hidiffm_%d%d',ti,tf);   % hori. diff.
end
if pn == 22 
filenameg = sprintf('wbudget_1_vidiffm_%d%d',ti,tf);   % vert. diff.
end
if pn == 23 
filenameg = sprintf('wbudget_1_nudiffm_%d%d',ti,tf);   % hori.+vert. numerical diff.
end
if pn == 31
filenameg = sprintf('wbudget_1_hturbm_%d%d',ti,tf);   % hori. turb.
end
if pn == 32
filenameg = sprintf('wbudget_1_vturbm_%d%d',ti,tf);   % vert. turb.
end
if pn == 4
filenameg = sprintf('wbudget_1_pgradm_%d%d',ti,tf);   % pgf
end
if pn == 5
filenameg = sprintf('wbudget_1_rdampm_%d%d',ti,tf);   % rdamp
end
if pn == 6
filenameg = sprintf('wbudget_1_buoym_%d%d',ti,tf);   % buoym.
end
if pn == 7
filenameg = sprintf('wbudget_1_myrhsm_%d%d',ti,tf);   % y-avg. r.h.s. terms
end
if pn == 8
filenameg = sprintf('wbudget_1_ddtwm_%d%d',ti,tf);   % y-avg. dwdt
end
if pn == 9
filenameg = sprintf('wbudget_1_profile_dwdt_rhs_%d%d',ti,tf);   % y-avg. dwdt
end
if pn == 10
filenameg = sprintf('wbudget_1_sgturb_%d%d',ti,tf);   % y-avg. dwdt
end

% saveas(f,fullfile(fnameg,filenameg),'epsc');













end   % end seq = 6:6

