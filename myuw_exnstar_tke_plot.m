close all
clear all

tic 

% directory for dwdt_eqn.mat files
% fname = '/aos/home/mtang/Documents/matrices';
% fname = 'D:\Users\stang33\Desktop';
fname = 'C:\Users\SiuLung\Downloads';   % for dx = 250

filename1 = 'myuw_exnstar_tke_dx125_may08.mat';
filename2 = 'myuw_exnstar_tke_dx250_mar07.mat';
filename5 = 'myuw_exnstar_tke_dx500_feb27.mat';

% directory for output graphs
% fnameg = '/aos/home/mtang/Documents/graphs';
% fnameg = 'D:\Users\stang33\Desktop';
fnameg = 'C:\Users\SiuLung\Downloads';

% input range of z interested
zmax = 4;
zmin = 0;

% input plot number interested
% (11) = y-avg. u, (12) = y-avg. w, 
% (2) = y-avg. exn-star
% (31) = subgrid tke, (32) = resolved tke, (33) = total tke
% (34) = profile of hori. mean (total) tke within region
% (4) = y-avg. pressure
% (20) = all

pn = 34;

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

% get y-avg. u matrices from arrays
myuhrm1 = terms1.terms.(fieldname).myuhrm;   % time avg. of y-avg. u for dx=125
myuhrm2 = terms2.terms.(fieldname).myuhrm;   % time avg. of y-avg. u for dx=250
myuhrm5 = terms5.terms.(fieldname).myuhrm;   % time avg. of y-avg. u for dx=500

% reduce myuhrm matrices into chosen region
myuhrm1rn = myuhrm1(:,zii:zfi);   % for dx=125
myuhrm2rn = myuhrm2(:,zii:zfi);   % for dx=250
myuhrm5rn = myuhrm5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of myuhrm1rn (to remove white space as z starts at 0.5)
myuhrm1r = [myuhrm1rn(:,1) myuhrm1rn];   % for dx=125
myuhrm2r = [myuhrm2rn(:,1) myuhrm2rn];   % for dx=250
myuhrm5r = [myuhrm5rn(:,1) myuhrm5rn];   % for dx=500

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

% get y-avg. exn-star matrices from arrays
myexnstarm1 = terms1.terms.(fieldname).myexnstarm;   % time avg. of y-avg. exn-star for dx=125
myexnstarm2 = terms2.terms.(fieldname).myexnstarm;   % time avg. of y-avg. exn-star for dx=250
myexnstarm5 = terms5.terms.(fieldname).myexnstarm;   % time avg. of y-avg. exn-star for dx=500

% reduce myexnstarm matrices into chosen region
myexnstarm1rn = myexnstarm1(:,zii:zfi);   % for dx=125
myexnstarm2rn = myexnstarm2(:,zii:zfi);   % for dx=250
myexnstarm5rn = myexnstarm5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of myexnstarm (to remove white space as z starts at 0.5)
myexnstarm1r = [myexnstarm1rn(:,1) myexnstarm1rn];   % for dx=125
myexnstarm2r = [myexnstarm2rn(:,1) myexnstarm2rn];   % for dx=250
myexnstarm5r = [myexnstarm5rn(:,1) myexnstarm5rn];   % for dx=500

% get y-avg. exn-star matrices from arrays
myexnstarm1 = terms1.terms.(fieldname).myexnstarm;   % time avg. of y-avg. exn-star for dx=125
myexnstarm2 = terms2.terms.(fieldname).myexnstarm;   % time avg. of y-avg. exn-star for dx=250
myexnstarm5 = terms5.terms.(fieldname).myexnstarm;   % time avg. of y-avg. exn-star for dx=500

% reduce myexnstarm matrices into chosen region
myexnstarm1rn = myexnstarm1(:,zii:zfi);   % for dx=125
myexnstarm2rn = myexnstarm2(:,zii:zfi);   % for dx=250
myexnstarm5rn = myexnstarm5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of myexnstarm (to remove white space as z starts at 0.5)
myexnstarm1r = [myexnstarm1rn(:,1) myexnstarm1rn];   % for dx=125
myexnstarm2r = [myexnstarm2rn(:,1) myexnstarm2rn];   % for dx=250
myexnstarm5r = [myexnstarm5rn(:,1) myexnstarm5rn];   % for dx=500

% get y-avg. rho*subgrid tke matrices from arrays
myrhosgtkem1 = terms1.terms.(fieldname).myrhosgtkem;   % time avg. of y-avg. rho*subgrid tke for dx=125
myrhosgtkem2 = terms2.terms.(fieldname).myrhosgtkem;   % time avg. of y-avg. rho*subgrid tke for dx=250
myrhosgtkem5 = terms5.terms.(fieldname).myrhosgtkem;   % time avg. of y-avg. rho*subgrid tke for dx=500

% reduce myrhosgtkem matrices into chosen region
myrhosgtkem1rn = myrhosgtkem1(:,zii:zfi);   % for dx=125
myrhosgtkem2rn = myrhosgtkem2(:,zii:zfi);   % for dx=250
myrhosgtkem5rn = myrhosgtkem5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of myrhosgtkem (to remove white space as z starts at 0.5)
myrhosgtkem1r = [myrhosgtkem1rn(:,1) myrhosgtkem1rn];   % for dx=125
myrhosgtkem2r = [myrhosgtkem2rn(:,1) myrhosgtkem2rn];   % for dx=250
myrhosgtkem5r = [myrhosgtkem5rn(:,1) myrhosgtkem5rn];   % for dx=500

% get y-avg. rho*resolved tke matrices from arrays
myrhoretkem1 = terms1.terms.(fieldname).myrhoretkem;   % time avg. of y-avg. rho*resolved tke for dx=125
myrhoretkem2 = terms2.terms.(fieldname).myrhoretkem;   % time avg. of y-avg. rho*resolved tke for dx=250
myrhoretkem5 = terms5.terms.(fieldname).myrhoretkem;   % time avg. of y-avg. rho*resolved tke for dx=500

% reduce myrhoretkem matrices into chosen region
myrhoretkem1rn = myrhoretkem1(:,zii:zfi);   % for dx=125
myrhoretkem2rn = myrhoretkem2(:,zii:zfi);   % for dx=250
myrhoretkem5rn = myrhoretkem5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of myrhoretkem (to remove white space as z starts at 0.5)
myrhoretkem1r = [myrhoretkem1rn(:,1) myrhoretkem1rn];   % for dx=125
myrhoretkem2r = [myrhoretkem2rn(:,1) myrhoretkem2rn];   % for dx=250
myrhoretkem5r = [myrhoretkem5rn(:,1) myrhoretkem5rn];   % for dx=500

% get y-avg. rho*total tke matrices from arrays
myrhottkem1 = terms1.terms.(fieldname).myrhottkem;   % time avg. of y-avg. rho*total tke for dx=125
myrhottkem2 = terms2.terms.(fieldname).myrhottkem;   % time avg. of y-avg. rho*total tke for dx=250
myrhottkem5 = terms5.terms.(fieldname).myrhottkem;   % time avg. of y-avg. rho*total tke for dx=500

% reduce myrhottkem matrices into chosen region
myrhottkem1rn = myrhottkem1(:,zii:zfi);   % for dx=125
myrhottkem2rn = myrhottkem2(:,zii:zfi);   % for dx=250
myrhottkem5rn = myrhottkem5(:,zii:zfi);   % for dx=500

% add one lvl. at the bottom of myrhottkem (to remove white space as z starts at 0.5)
myrhottkem1r = [myrhottkem1rn(:,1) myrhottkem1rn];   % for dx=125
myrhottkem2r = [myrhottkem2rn(:,1) myrhottkem2rn];   % for dx=250
myrhottkem5r = [myrhottkem5rn(:,1) myrhottkem5rn];   % for dx=500


%--------------------------------------------------------------------------------
% set axes font size
% axfs = 12;
tfs = 14;   % title font size
tefs = 14;   % text font size (plot number)
lbfs = 13;   % x,y-label font size

%-------------------------
% original range of colorbar
% if pn == 11 || pn == 12 
% v250 = -5:5/10:5;   % for dx=250
% vdiff = -0.6:0.6/6:0.6;   % for dx=125-dx=250 and dx=500-dx=250
% end
% 
% if pn == 2
% v250 = -6e-5:6e-5/10:6e-5;   % for dx=250 
% vdiff = -6e-6:6e-6/6:6e-6;
% end
% 
% if pn == 31 || pn == 32 || pn == 33
% v250 = -2.5:2.5/10:2.5;   % for dx=250
% vdiff = -1.5:1.5/10:1.5;   % for dx=125-dx=250 and dx=500-dx=250
% end
% % 
% if pn == 31
% v250 = -0.8:0.8/10:0.8;
% vdiff = -0.2:0.2/10:0.2;
% end
%-------------------------
%-------------------------
% change range of colorbar to 0.5* that of mean plots

if pn == 11 || pn == 12 
v250 = -5:5/10:5;   % for dx=250
vdiff = 0.5.*v250;   % for dx=125-dx=250 and dx=500-dx=250
end

if pn == 2
v250 = -6e-5:6e-5/10:6e-5;   % for dx=250 
vdiff = 0.5.*v250;
end

% if pn == 31 || pn == 32 || pn == 33
% v250 = -5:5/10:5;   % for dx=250
% vdiff = 0.5.*v250;   % for dx=125-dx=250 and dx=500-dx=250
% end

if pn == 31 || pn == 32 || pn == 33
% v250 = -3.15:0.1:3.15;
% vdiff = 0.5.*v250;   % for dx=125-dx=250 and dx=500-dx=250

% v250 = -3.15:0.1:3.15;
v250 = -3.225:0.15:3.225;
vdiff = 0.5.*v250;   % for dx=125-dx=250 and dx=500-dx=250
end

llvp = 0.25:0.5:3.75;   % for positive y-avg. w contours
llvn = -1.75:0.5:-0.25;   % for negative y-avg. w contours   % for negative y-avg. w contours

dllvp = 0.1:0.2:0.5;   % for diff. in positive y-avg. w contours
dllvn = -0.9:0.2:-0.1;   % for diff. in negative y-avg. w contours


xmina = -3;   % xmin for all plots
xmaxa = 3;   % xmax for all plots
xticka = 1.5;   % xticks for all plots

pbra = 2;   % apsect ratio for plots, y:x ratio

%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------

% (I_A) MAKE y-avg. u PLOT

if pn == 11 || pn == 20
     
% create vector definiing contour lines
vmyu = v250;

% figure
f = figure('units','normalized','outerposition',[0 0 1 1]);

%---------------------------------------
% for dx=250, place it in the middle

% To prevent white patches formed on contour plot because of entries in
% myuhrm2r exceeding limits set by v250,
% find those entries and set them as the limit of v250

% initialize new myuhrm2r
myuhrm2rnew = myuhrm2r;

if find(myuhrm2r > v250(end)) ~= 0   % for upper limit of v250
    myuhrm2rnew(find(myuhrm2r > v250(end))) = v250(end);
end

if find(myuhrm2r < v250(1)) ~= 0   % for lower limit of v250
    myuhrm2rnew(find(myuhrm2r < v250(1))) = v250(1);
end



subplot(1,3,2);
% [C, h1] = contourf(xh2d2,z2d2,myuhrm2rnew,vmyu);
[C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),myuhrm2rnew(220:260,:),vmyu);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyu(1) vmyu(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\overline{u}\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
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

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
% for dx=125 - dx=250

% interpolate myuhrm from dx=125 to dx=250 case
myuhrm1r2 = 0.5.*(myuhrm1r(1:2:end-1,:)+ myuhrm1r(2:2:end,:));

% get difference btw. interpolated myuhrm1r2 from myuhrm2r
diffmyuhrm12 = myuhrm1r2 - myuhrm2r;

% set contour lines 
vmyu12 = vdiff;

% To prevent white patches formed on contour plot because of entries in
% diffmyuhrm12 exceeding limits set by vdiff,
% find those entries and set them as the limit of vdiff

% initialize new diffmyuhrm12
diffmyuhrm12new = diffmyuhrm12;

if find(diffmyuhrm12 > vdiff(end)) ~= 0   % for upper limit of vdiff
    diffmyuhrm12new(find(diffmyuhrm12 > vdiff(end))) = vdiff(end);
end

if find(diffmyuhrm12 < vdiff(1)) ~= 0   % for lower limit of vdiff
    diffmyuhrm12new(find(diffmyuhrm12 < vdiff(1))) = 0.5*(vdiff(1)+vdiff(2));
end

subplot(1,3,1);
[C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),diffmyuhrm12new(220:260,:),vmyu12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyu12(1) vmyu12(end)];   % range of caxis
caxis([cr(1) cr(end)]);

cbh = colorbar;   % create colobar
cbh.Ticks = vdiff(1:2:end);   % set ticks on colobar
cbh.TickLabels = vdiff(1:2:end);   % set tick labels on colobar

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

text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 

% end for dx=125
%---------------------------------------
% for dx=500 - dx=250

% interpolate myuhrm from dx=250 to dx=500 case
myuhrm2r5 = 0.5.*(myuhrm2r(1:2:end-1,:)+ myuhrm2r(2:2:end,:));

% get difference btw. interpolated myuhrm5r from myuhrm2r5
diffmyuhrm52 = myuhrm5r - myuhrm2r5;

% set contour lines 
vmyu12 = vdiff;

% To prevent white patches formed on contour plot because of entries in
% diffmyuhrm52 exceeding limits set by vdiff,
% find those entries and set them as the limit of vdiff

% initialize new diffmyuhrm12
diffmyuhrm52new = diffmyuhrm52;

if find(diffmyuhrm52 > vdiff(end)) ~= 0   % for upper limit of vdiff
    diffmyuhrm52new(find(diffmyuhrm52 > vdiff(end))) = vdiff(end);
end

if find(diffmyuhrm52 < vdiff(1)) ~= 0   % for lower limit of vdiff
    diffmyuhrm52new(find(diffmyuhrm52 < vdiff(1))) = 0.5*(vdiff(1)+vdiff(2));
end

subplot(1,3,3);
[C, h1] = contourf(xh2d5(100:140,:),z2d5(100:140,:),diffmyuhrm52new(100:140,:),vmyu12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyu12(1) vmyu12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbh = colorbar;   % create colobar
cbh.Ticks = vdiff(1:2:end);   % set ticks on colobar
cbh.TickLabels = vdiff(1:2:end);   % set tick labels on colobar

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

text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 

end   % end if pn == 91

% end for dx=500
%---------------------------------------
% END MAKING y-avg. u PLOT
%--------------------------------------------------------------------------------

% (I_B) MAKE y-avg. w PLOT

if pn == 12 || pn == 20
     
% create vector definiing contour lines
vmyw = v250;

% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle

% To prevent white patches formed on contour plot because of entries in
% myuhrm2r exceeding limits set by v250,
% find those entries and set them as the limit of v250

% initialize new myuhrm2r
mywhrm2rnew = mywhrm2r;

if find(mywhrm2r > v250(end)) ~= 0   % for upper limit of v250
    mywhrm2rnew(find(mywhrm2r > v250(end))) = v250(end);
end

if find(mywhrm2r < v250(1)) ~= 0   % for lower limit of v250
    mywhrm2rnew(find(mywhrm2r < v250(1))) = v250(1);
end


subplot(1,3,2);
[C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),mywhrm2rnew(220:260,:),vmyw);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyw(1) vmyw(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\overline{w}\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
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

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
% for dx=125 - dx=250

% interpolate mywhrm from dx=125 to dx=250 case
mywhrm1r2 = 0.5.*(mywhrm1r(1:2:end-1,:)+ mywhrm1r(2:2:end,:));

% get difference btw. interpolated mywhrm1r2 from mywhrm2r
diffmywhrm12 = mywhrm1r2 - mywhrm2r;

% set contour lines 
vmyw12 = vdiff;

% To prevent white patches formed on contour plot because of entries in
% diffmyuhrm12 exceeding limits set by vdiff,
% find those entries and set them as the limit of vdiff

% initialize new diffmywhrm12
diffmywhrm12new = diffmywhrm12;

if find(diffmywhrm12 > vdiff(end)) ~= 0   % for upper limit of vdiff
    diffmywhrm12new(find(diffmywhrm12 > vdiff(end))) = vdiff(end);
end

if find(diffmywhrm12 < vdiff(1)) ~= 0   % for lower limit of vdiff
    diffmywhrm12new(find(diffmywhrm12 < vdiff(1))) = 0.5*(vdiff(1) + vdiff(2));
end

subplot(1,3,1);
[C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),diffmywhrm12new(220:260,:),vmyw12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyw12(1) vmyw12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbh = colorbar;   % create colobar
cbh.Ticks = vdiff(1:2:end);   % set ticks on colobar
cbh.TickLabels = vdiff(1:2:end);   % set tick labels on colobar

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

text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 

% end for dx=125
%---------------------------------------
% for dx=500 - dx=250

% interpolate mywhrm from dx=250 to dx=500 case
mywhrm2r5 = 0.5.*(mywhrm2r(1:2:end-1,:)+ mywhrm2r(2:2:end,:));

% get difference btw. interpolated mywhrm5r from myuhrm2r5
diffmywhrm52 = mywhrm5r - mywhrm2r5;

% set contour lines 
vmyw12 = vdiff;

subplot(1,3,3);
[C, h1] = contourf(xh2d5(100:140,:),z2d5(100:140,:),diffmywhrm52(100:140,:),vmyw12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyw12(1) vmyw12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbh = colorbar;   % create colobar
cbh.Ticks = vdiff(1:2:end);   % set ticks on colobar
cbh.TickLabels = vdiff(1:2:end);   % set tick labels on colobar

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

text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 

end   % end if pn == 91

% end for dx=500
%---------------------------------------
% END MAKING y-avg. w PLOT
%--------------------------------------------------------------------------------

% (II) MAKE y-avg. exn-star PLOT

if pn == 2 || pn == 20
     
% create vector definiing contour lines
vmyes = v250;

% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle

% % To prevent white patches formed on contour plot because of entries in
% % myexnstarm2r exceeding limits set by v250,
% % find those entries and set them as the limit of v250
% 
% % initialize new myexnstarm2r
% myexnstarm2rnew = myexnstarm2r;
% 
% if find(myexnstarm2r > v250(end)) ~= 0   % for upper limit of v250
%     myexnstarm2rnew(find(myexnstarm2r > v250(end))) = v250(end);
% end
% 
% if find(myexnstarm2r < v250(1)) ~= 0   % for lower limit of v250
%     myexnstarm2rnew(find(myexnstarm2r < v250(1))) = v250(1);
% end

% To remove sheet-like appearance in contour plot of y-avg. exn-star at dx=250

% get horizontal avg. of myexnstarm2r
hmmyexnstarm2r = zeros(1,size(myexnstarm2r,2));
for j = 1:size(myexnstarm2r,2)
    hmmyexnstarm2r(j) = mean(myexnstarm2r(:,j));
end

% repmat hori. avg. of myexnstarm2r into 2-d matrix
hmmyexnstarm2r2d = repmat(hmmyexnstarm2r,[size(myexnstarm2r,1) 1]);

% subtract hori. avg. myexnstarm2r from myexnstarm2r
myexnstarm2hm = myexnstarm2r - hmmyexnstarm2r2d;


subplot(1,3,2);
% [C, h1] = contourf(xh2d2,z2d2,myexnstarm2r,vmyes);
[C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),myexnstarm2hm(220:260,:),vmyes);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyes(1) vmyes(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\overline{\pi^{*}}\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
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

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
% for dx=125 - dx=250

% interpolate myexnstarm1r from dx=125 to dx=250 case
myexnstarm1r2 = 0.5.*(myexnstarm1r(1:2:end-1,:)+ myexnstarm1r(2:2:end,:));

% get difference btw. interpolated myexnstarm1r2 from mywhrm2r
diffmyexnstarm12 = myexnstarm1r2 - myexnstarm2r;

% set contour lines 
vmyes12 = vdiff;

% To remove sheet-like appearance in contour plot of y-avg. exn-star at
% dx=125-dx=250

% get horizontal avg. of diffmyexnstarm12
hmdiffmyexnstarm12 = zeros(1,size(diffmyexnstarm12,2));
for j = 1:size(diffmyexnstarm12,2)
    hmdiffmyexnstarm12(j) = mean(diffmyexnstarm12(:,j));
end

% repmat hori. avg. of diffmyexnstarm12 into 2-d matrix
hmdiffmyexnstarm122d = repmat(hmdiffmyexnstarm12,[size(diffmyexnstarm12,1) 1]);

% subtract hori. avg. diffmyexnstarm12 from hmdiffmyexnstarm122d
diffmyexnstarm12hm = diffmyexnstarm12 - hmdiffmyexnstarm122d;

% To prevent white patches formed on contour plot because of entries in
% diffmyexnstarm12hm exceeding limits set by vdiff,
% find those entries and set them as the limit of vdiff

% initialize new diffmyexnstarm12hm
diffmyexnstarm12hmnew = diffmyexnstarm12hm;

if find(diffmyexnstarm12hm > vdiff(end)) ~= 0   % for upper limit of vdiff
    diffmyexnstarm12hmnew(find(diffmyexnstarm12hm > vdiff(end))) = vdiff(end);
end

if find(diffmyexnstarm12hm < vdiff(1)) ~= 0   % for lower limit of vdiff
    diffmyexnstarm12hmnew(find(diffmyexnstarm12hm < vdiff(1))) = 0.5*(vdiff(1)+vdiff(2));
end


subplot(1,3,1);
[C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),diffmyexnstarm12hmnew(220:260,:),vmyes12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyes12(1) vmyes12(end)];   % range of caxis
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

text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 

% end for dx=125
%---------------------------------------
% for dx=500 - dx=250

% interpolate myexnstarm2r from dx=250 to dx=500 case
myexnstarm2r5 = 0.5.*(myexnstarm2r(1:2:end-1,:)+ myexnstarm2r(2:2:end,:));

% get difference btw. interpolated mywhrm5r from myuhrm2r5
diffmyexnstarm52 = myexnstarm5r - myexnstarm2r5;

% set contour lines 
vmyes12 = vdiff;

% To remove sheet-like appearance in contour plot of y-avg. exn-star at
% dx=500-dx=250

% get horizontal avg. of diffmyexnstarm52
hmdiffmyexnstarm52 = zeros(1,size(diffmyexnstarm52,2));
for j = 1:size(diffmyexnstarm52,2)
    hmdiffmyexnstarm52(j) = mean(diffmyexnstarm52(:,j));
end

% repmat hori. avg. of diffmyexnstarm52 into 2-d matrix
hmdiffmyexnstarm522d = repmat(hmdiffmyexnstarm52,[size(diffmyexnstarm52,1) 1]);

% subtract hori. avg. diffmyexnstarm52 from hmdiffmyexnstarm522d
diffmyexnstarm52hm = diffmyexnstarm52 - hmdiffmyexnstarm522d;

subplot(1,3,3);
[C, h1] = contourf(xh2d5(100:140,:),z2d5(100:140,:),diffmyexnstarm52hm(100:140,:),vmyes12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyes12(1) vmyes12(end)];   % range of caxis
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

text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 

end   % end if pn == 2

% end for dx=500
%---------------------------------------
% END MAKING y-avg. exn-star PLOT
%--------------------------------------------------------------------------------

% (III_A) MAKE y-avg. rho* subgrid tke PLOT

if pn == 31 || pn == 20
     
% create vector definiing contour lines
vmyrsgtke = v250;

% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle

% % To prevent white patches formed on contour plot because of entries in
% % myexnstarm2r exceeding limits set by v250,
% % find those entries and set them as the limit of v250
% 
% % initialize new myexnstarm2r
% myexnstarm2rnew = myexnstarm2r;
% 
% if find(myexnstarm2r > v250(end)) ~= 0   % for upper limit of v250
%     myexnstarm2rnew(find(myexnstarm2r > v250(end))) = v250(end);
% end
% 
% if find(myexnstarm2r < v250(1)) ~= 0   % for lower limit of v250
%     myexnstarm2rnew(find(myexnstarm2r < v250(1))) = v250(1);
% end

% To remove sheet-like appearance in contour plot of y-avg. exn-star at dx=250

% get horizontal avg. of myrhosgtkem2r
hmmyrhosgtkem2r = zeros(1,size(myrhosgtkem2r,2));
for j = 1:size(myrhosgtkem2r,2)
    hmmyrhosgtkem2r(j) = mean(myrhosgtkem2r(:,j));
end

% repmat hori. avg. of myrhosgtkem2r into 2-d matrix
hmmyrhosgtkem2r2d = repmat(hmmyrhosgtkem2r,[size(myrhosgtkem2r,1) 1]);

% subtract hori. avg. myrhosgtkem2r from myrhosgtkem2r
myrhosgtkem2hm = myrhosgtkem2r - hmmyrhosgtkem2r2d;


subplot(1,3,2);
[C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),myrhosgtkem2r(220:260,:),vmyrsgtke);
% [C, h1] = contourf(xh2d2,z2d2,myrhosgtkem2hm,vmyrsgtke);


cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyrsgtke(1) vmyrsgtke(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\overline{\textrm{subgrid tke}}\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
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

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
% for dx=125 - dx=250

% interpolate myexnstarm1r from dx=125 to dx=250 case
myrhosgtkem1r2 = 0.5.*(myrhosgtkem1r(1:2:end-1,:)+ myrhosgtkem1r(2:2:end,:));

% get difference btw. interpolated myrhosgtkem1r2 from myrhosgtkem2r
diffmyrhosgtkem12 = myrhosgtkem1r2 - myrhosgtkem2r;

% set contour lines 
vmyrsgtke12 = vdiff;

% % To remove sheet-like appearance in contour plot of y-avg. exn-star at
% % dx=125-dx=250
% 
% % get horizontal avg. of diffmyexnstarm12
% hmdiffmyexnstarm12 = zeros(1,size(diffmyexnstarm12,2));
% for j = 1:size(diffmyexnstarm12,2)
%     hmdiffmyexnstarm12(j) = mean(diffmyexnstarm12(:,j));
% end
% 
% % repmat hori. avg. of diffmyexnstarm12 into 2-d matrix
% hmdiffmyexnstarm122d = repmat(hmdiffmyexnstarm12,[size(diffmyexnstarm12,1) 1]);
% 
% % subtract hori. avg. diffmyexnstarm12 from hmdiffmyexnstarm122d
% diffmyexnstarm12hm = diffmyexnstarm12 - hmdiffmyexnstarm122d;

% To prevent white patches formed on contour plot because of entries in
% diffmyrhosgtkem12 exceeding limits set by vdiff,
% find those entries and set them as the limit of vdiff

% initialize new diffmyrhosgtkem12
diffmyrhosgtkem12new = diffmyrhosgtkem12;

if find(diffmyrhosgtkem12 > vdiff(end)) ~= 0   % for upper limit of v250
    diffmyrhosgtkem12new(find(diffmyrhosgtkem12 > vdiff(end))) = vdiff(end);
end

if find(diffmyrhosgtkem12 < vdiff(1)) ~= 0   % for lower limit of v250
    diffmyrhosgtkem12new(find(diffmyrhosgtkem12 < vdiff(1))) = 0.5*(vdiff(1) + vdiff(2));
end

subplot(1,3,1);
[C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),diffmyrhosgtkem12new(220:260,:),vmyrsgtke12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyrsgtke12(1) vmyrsgtke12(end)];   % range of caxis
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

text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 

% end for dx=125
%---------------------------------------
% for dx=500 - dx=250

% interpolate myrhosgtkem2r from dx=250 to dx=500 case
myrhosgtkem2r5 = 0.5.*(myrhosgtkem2r(1:2:end-1,:)+ myrhosgtkem2r(2:2:end,:));

% get difference btw. interpolated myrhosgtkem5r from myrhosgtkem2r5
diffmyrhosgtkem52 = myrhosgtkem5r - myrhosgtkem2r5;

% set contour lines 
vmyrsgtke12 = vdiff;

% To remove sheet-like appearance in contour plot of y-avg. exn-star at
% dx=500-dx=250

% get horizontal avg. of diffmyrhosgtkem52
hmdiffmyrhosgtkem52 = zeros(1,size(diffmyrhosgtkem52,2));
for j = 1:size(diffmyrhosgtkem52,2)
    hmdiffmyrhosgtkem52(j) = mean(diffmyrhosgtkem52(:,j));
end

% repmat hori. avg. of diffmyrhosgtkem52 into 2-d matrix
hmdiffmyrhosgtkem522d = repmat(hmdiffmyrhosgtkem52,[size(diffmyrhosgtkem52,1) 1]);

% subtract hori. avg. diffmyrhosgtkem52 from hmdiffmyrhosgtkem522d
diffmyrhosgtkem52hm = diffmyrhosgtkem52 - hmdiffmyrhosgtkem522d;

subplot(1,3,3);
[C, h1] = contourf(xh2d5(100:140,:),z2d5(100:140,:),diffmyrhosgtkem52(100:140,:),vmyrsgtke12);
% [C, h1] = contourf(xh2d5,z2d5,diffmyrhosgtkem52hm,vmyrsgtke12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyrsgtke12(1) vmyrsgtke12(end)];   % range of caxis
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

text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 

end   % end if pn == 2

% end for dx=500
%---------------------------------------
% END MAKING y-avg. subgrid tke PLOT
%--------------------------------------------------------------------------------

% (III_B) MAKE y-avg. rho* resolved tke PLOT

if pn == 32 || pn == 20
     
% create vector definiing contour lines
vmyrretke = v250;

% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle

% % To prevent white patches formed on contour plot because of entries in
% % myexnstarm2r exceeding limits set by v250,
% % find those entries and set them as the limit of v250
% 
% % initialize new myexnstarm2r
% myexnstarm2rnew = myexnstarm2r;
% 
% if find(myexnstarm2r > v250(end)) ~= 0   % for upper limit of v250
%     myexnstarm2rnew(find(myexnstarm2r > v250(end))) = v250(end);
% end
% 
% if find(myexnstarm2r < v250(1)) ~= 0   % for lower limit of v250
%     myexnstarm2rnew(find(myexnstarm2r < v250(1))) = v250(1);
% end

% To remove sheet-like appearance in contour plot of y-avg. exn-star at dx=250

% get horizontal avg. of myrhoretkem2r
hmmyrhoretkem2r = zeros(1,size(myrhoretkem2r,2));
for j = 1:size(myrhoretkem2r,2)
    hmmyrhoretkem2r(j) = mean(myrhoretkem2r(:,j));
end

% repmat hori. avg. of myrhoretkem2r into 2-d matrix
hmmyrhoretkem2r2d = repmat(hmmyrhoretkem2r,[size(myrhoretkem2r,1) 1]);

% subtract hori. avg. myrhosgtkem2r from myrhosgtkem2r
myrhoretkem2hm = myrhoretkem2r - hmmyrhoretkem2r2d;


subplot(1,3,2);
[C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),myrhoretkem2r(220:260,:),vmyrretke);
% [C, h1] = contourf(xh2d2,z2d2,myrhoretkem2hm,vmyrretke);


cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyrretke(1) vmyrretke(end)];   % range of caxis
caxis([cr(1) cr(end)]);
colorbar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\overline{\textrm{resolved tke}}\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
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

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% end for dx=250
%---------------------------------------
% for dx=125 - dx=250

% interpolate myrhosgtkem1r from dx=125 to dx=250 case
myrhoretkem1r2 = 0.5.*(myrhoretkem1r(1:2:end-1,:)+ myrhoretkem1r(2:2:end,:));

% get difference btw. interpolated myrhoretkem1r2 from myrhoretkem2r
diffmyrhoretkem12 = myrhoretkem1r2 - myrhoretkem2r;

% set contour lines 
vmyrretke12 = vdiff;

% To remove sheet-like appearance in contour plot of y-avg. exn-star at
% dx=125-dx=250

% get horizontal avg. of diffmyexnstarm12
hmdiffmyrhoretkem12 = zeros(1,size(diffmyrhoretkem12,2));
for j = 1:size(diffmyrhoretkem12,2)
    hmdiffmyrhoretkem12(j) = mean(diffmyrhoretkem12(:,j));
end

% repmat hori. avg. of diffmyrhoretkem12 into 2-d matrix
hmdiffmyrhoretkem122d = repmat(hmdiffmyrhoretkem12,[size(diffmyrhoretkem12,1) 1]);

% subtract hori. avg. diffmyrhoretkem12 from hmdiffmyrhoretkem122d
diffmyrhoretkem12hm = diffmyrhoretkem12 - hmdiffmyrhoretkem122d;


subplot(1,3,1);
[C, h1] = contourf(xh2d2(220:260,:),z2d2(220:260,:),diffmyrhoretkem12(220:260,:),vmyrretke12);
% [C, h1] = contourf(xh2d2,z2d2,diffmyrhoretkem12hm,vmyrretke12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyrretke12(1) vmyrretke12(end)];   % range of caxis
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

text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 

% end for dx=125
%---------------------------------------
% for dx=500 - dx=250

% interpolate myrhoretkem2r from dx=250 to dx=500 case
myrhoretkem2r5 = 0.5.*(myrhoretkem2r(1:2:end-1,:)+ myrhoretkem2r(2:2:end,:));

% get difference btw. interpolated myrhoretkem5r from myrhoretkem2r5
diffmyrhoretkem52 = myrhoretkem5r - myrhoretkem2r5;

% set contour lines 
vmyrretke12 = vdiff;

% To remove sheet-like appearance in contour plot of y-avg. exn-star at
% dx=500-dx=250

% get horizontal avg. of diffmyrhoretkem52
hmdiffmyrhoretkem52 = zeros(1,size(diffmyrhoretkem52,2));
for j = 1:size(diffmyrhoretkem52,2)
    hmdiffmyrhoretkem52(j) = mean(diffmyrhoretkem52(:,j));
end

% repmat hori. avg. of diffmyrhoretkem52 into 2-d matrix
hmdiffmyrhoretkem522d = repmat(hmdiffmyrhoretkem52,[size(diffmyrhoretkem52,1) 1]);

% subtract hori. avg. diffmyrhoretkem52 from hmdiffmyrhoretkem522d
diffmyrhoretkem52hm = diffmyrhoretkem52 - hmdiffmyrhoretkem522d;

subplot(1,3,3);
[C, h1] = contourf(xh2d5(100:140,:),z2d5(100:140,:),diffmyrhoretkem52(100:140,:),vmyrretke12);
% [C, h1] = contourf(xh2d5,z2d5,diffmyrhoretkem52hm,vmyrretke12);

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyrretke12(1) vmyrretke12(end)];   % range of caxis
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

text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 

end   % end if pn == 2

% end for dx=500
%---------------------------------------
% END MAKING y-avg. resolved tke PLOT
%--------------------------------------------------------------------------------

% (III_C) MAKE y-avg. rho* total tke PLOT

if pn == 33 || pn == 20
     
% create vector definiing contour lines
vmyrttke = v250;

% figure
f = figure('units','normalized','outerposition',[0 0 1 1])

%---------------------------------------
% for dx=250, place it in the middle

% % To prevent white patches formed on contour plot because of entries in
% % myexnstarm2r exceeding limits set by v250,
% % find those entries and set them as the limit of v250
% 
% % initialize new myexnstarm2r
% myexnstarm2rnew = myexnstarm2r;
% 
% if find(myexnstarm2r > v250(end)) ~= 0   % for upper limit of v250
%     myexnstarm2rnew(find(myexnstarm2r > v250(end))) = v250(end);
% end
% 
% if find(myexnstarm2r < v250(1)) ~= 0   % for lower limit of v250
%     myexnstarm2rnew(find(myexnstarm2r < v250(1))) = v250(1);
% end

% To remove sheet-like appearance in contour plot of y-avg. exn-star at dx=250

% get horizontal avg. of myrhottkem2r
hmmyrhottkem2r = zeros(1,size(myrhottkem2r,2));
for j = 1:size(myrhottkem2r,2)
    hmmyrhottkem2r(j) = mean(myrhottkem2r(:,j));
end

% repmat hori. avg. of myrhottkem2r into 2-d matrix
hmmyrhottkem2r2d = repmat(hmmyrhottkem2r,[size(myrhottkem2r,1) 1]);

% subtract hori. avg. myrhottkem2r from myrhottkem2r
myrhottkem2hm = myrhottkem2r - hmmyrhottkem2r2d;


subplot(3,1,2);
[C12, h12] = contourf(xh2d2(220:260,1:35),z2d2(220:260,1:35),myrhottkem2r(220:260,1:35),vmyrttke); hold on
[C22, h22] = contour(xh2d2(220:260,1:35),z2d2(220:260,1:35),mywhrm2r(220:260,1:35),llvp); hold on
[C32, h32] = contour(xh2d2(220:260,1:35),z2d2(220:260,1:35),mywhrm2r(220:260,1:35),llvn); hold off

% set line style and color on each contour object of y-avg. w
h22.LineStyle = '-';   % solid line for positive y-avg. w contours
h22.LineColor = 'k';   % set contour lines at black

h32.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h32.LineColor = 'k';   % set contour lines at black

% set line color in contours of y-avg. total tke
h12.LineColor = 'none';   % set dudx contour lines as colorless

grid on 
box on

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyrttke(1) vmyrttke(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbh2 = colorbar;   % create colorbar
cbh2.Position(1) = 0.9;
cbh2.Ticks = -3:1:3;   % set ticks on colobar
% cbh2.Ticks = v250(1):0.5:v250(end);   % set ticks on colobar

 
% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(\overline{\textrm{total tke}}\), ' num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

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

% pbr = pbra;
pbr = 1/3;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

text(-2.85,2.6,'(b)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % upper left corner
text(2.05,2.7,'\(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=250
%---------------------------------------
% for dx=125 - dx=250

%-----------------------
% use interp1 to interpolate grid points from 250-m to 125-m case
myrhottkem2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));   % y-avg. hori. subgrid turb. flux of w' from dx=250 to dx=125
mywhrm2r1 = zeros(size(mywhrm1r,1),size(mywhrm1r,2));  % y-avg. w from dx=250 to dx=125

for k = 1:size(mywhrm1r,2)
myrhottkem2r1(:,k) = interp1(xh2(xhii2:xhfi2),myrhottkem2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
mywhrm2r1(:,k) = interp1(xh2(xhii2:xhfi2),mywhrm2r(:,k),xh1(xhii1:xhfi1),'linear','extrap');
end

% get diff. between 125-m and 250-m cases
diffmyrhottkem12 = myrhottkem1r - myrhottkem2r1;
diffmywhrm12 = mywhrm1r - mywhrm2r1;   % diff. in w (dx=125-dx=250)

%-----------------------

% set contour lines 
vmyrttke12 = vdiff;

% To prevent white space formed because value at that point exceeds
% colorbar limits

% initialize new diffmyrhottkem12
diffmyrhottkem12new = diffmyrhottkem12;

if find(diffmyrhottkem12> vmyrttke12(end)) ~= 0   % for upper limit of v250
    diffmyrhottkem12new(find(diffmyrhottkem12 > vmyrttke12(end))) = vmyrttke12(end);
end

if find(diffmyrhottkem12 < vmyrttke12(1)) ~= 0   % for lower limit of v250
    diffmyrhottkem12new(find(diffmyrhottkem12 < vmyrttke12(1))) = vmyrttke12(1);
end


subplot(3,1,1);
[C11, h11] = contourf(xh2d1(450:510,1:35),z2d1(450:510,1:35),diffmyrhottkem12new(450:510,1:35),vmyrttke12); hold on
[C21, h21] = contour(xh2d1(450:510,1:35),z2d1(450:510,1:35),diffmywhrm12(450:510,1:35),dllvp); hold on
[C31, h31] = contour(xh2d1(450:510,1:35),z2d1(450:510,1:35),diffmywhrm12(450:510,1:35),dllvn); hold off

% set line style and color on each contour object of y-avg. w
h21.LineStyle = '-';   % solid line for positive y-avg. w contours
h21.LineColor = 'k';   % set contour lines at black

h31.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h31.LineColor = 'k';   % set contour lines at black

% set line color in contours of y-avg. total tke
h11.LineColor = 'none';   % set dudx contour lines as colorless

grid on 
box on

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyrttke12(1) vmyrttke12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbh1 = colorbar;   % create colorbar
cbh1.Position(1) = 0.9;
cbh1.Ticks = -1.5:0.5:1.5;   % set ticks on colobar
% cbh1.Ticks = vdiff(1):0.5:vdiff(end);   % set ticks on colobar

% xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

title(['\(\overline{\textrm{Total TKE}}\)'],'Interpreter','latex','FontSize',tfs);

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

% pbr = pbra;
pbr = 1/3;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(-2.85,2.6,'(a)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % upper left corner
text(0.9,2.7,'\(\Delta=125\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end for dx=125
%---------------------------------------
% for dx=500 - dx=250

%-----------------------
% use interp1 to interpolate grid points from 500-m to 250-m case
myrhottkem5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));   % y-avg. hori. subgrid turb. flux of w' from dx=500 to dx=250
mywhrm5r2 = zeros(size(mywhrm2r,1),size(mywhrm2r,2));  % y-avg. w from dx=500 to dx=250

for k = 1:size(mywhrm2r,2)
myrhottkem5r2(:,k) = interp1(xh5(xhii5:xhfi5),myrhottkem5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
mywhrm5r2(:,k) = interp1(xh5(xhii5:xhfi5),mywhrm5r(:,k),xh2(xhii2:xhfi2),'linear','extrap');
end

% get diff. between 500-m and 250-m cases
diffmyrhottkem52 = myrhottkem5r2 - myrhottkem2r;
diffmywhrm52 = mywhrm5r2 - mywhrm2r;   % diff. in w (dx=500-dx=250)

%-----------------------

% set contour lines 
vmyrttke12 = vdiff;

% To remove sheet-like appearance in contour plot of y-avg. exn-star at
% dx=500-dx=250

% get horizontal avg. of diffmyrhottkem52
hmdiffmyrhottkem52 = zeros(1,size(diffmyrhottkem52,2));
for j = 1:size(diffmyrhottkem52,2)
    hmdiffmyrhottkem52(j) = mean(diffmyrhottkem52(:,j));
end

% repmat hori. avg. of diffmyrhottkem52 into 2-d matrix
hmdiffmyrhottkem522d = repmat(hmdiffmyrhottkem52,[size(diffmyrhottkem52,1) 1]);

% subtract hori. avg. diffmyrhottkem52 from hmdiffmyrhottkem522d
diffmyrhottkem52hm = diffmyrhottkem52 - hmdiffmyrhottkem522d;

subplot(3,1,3);
[C15, h15] = contourf(xh2d2(220:260,1:35),z2d2(220:260,1:35),diffmyrhottkem52(220:260,1:35),vmyrttke12); hold on
[C25, h25] = contour(xh2d2(220:260,1:35),z2d2(220:260,1:35),diffmywhrm52(220:260,1:35),dllvp); hold on
[C35, h35] = contour(xh2d2(220:260,1:35),z2d2(220:260,1:35),diffmywhrm52(220:260,1:35),dllvn); hold off

% set line style and color on each contour object of y-avg. w
h25.LineStyle = '-';   % solid line for positive y-avg. w contours
h25.LineColor = 'k';   % set contour lines at black

h35.LineStyle = '--';   % dashed line for negative y-avg. w contours 
h35.LineColor = 'k';   % set contour lines at black

% set line color in contours of y-avg. total tke
h15.LineColor = 'none';   % set dudx contour lines as colorless

grid on 
box on

cm = load('cmap_bous.mat');   %load edited color map
cmap = cm.cmap;          %access variable from structure
colormap(cmap);   %insert chosen colormap

cr = [vmyrttke12(1) vmyrttke12(end)];   % range of caxis
caxis([cr(1) cr(end)]);
cbh5 = colorbar
cbh5.Position(1) = 0.9;
cbh5.Ticks = -1.5:0.5:1.5;   % set ticks on colobar

xlabel('\textit{x}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',lbfs);

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

% pbr = pbra;
pbr = 1/3;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(3.4,3.16,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold'); 
text(-2.85,2.6,'(c)', 'Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');   % upper left corner
text(0.9,2.7,'\(\Delta=500\) m \(-\) \(\Delta=250\) m', 'Interpreter', 'latex','FontSize',tefs-4,'FontWeight','bold');

% end   % end if pn == 33

% end for dx=500
%---------------------------------------

% SET SUBPLOT POSITIONS

% reduce vert. spacing between subplots
ax2.Position(2) = ax5.Position(2) + ax5.Position(4) + 0.0425;
ax1.Position(2) = ax2.Position(2) + ax2.Position(4) + 0.0425;

% set colorbar width
cbh1.Position(3) = 0.5*cbh1.Position(3);
cbh2.Position(3) = 0.5*cbh2.Position(3);
cbh5.Position(3) = 0.5*cbh5.Position(3);

% reduce hori. spacing between subplots and colorbars
cbh1.Position(1) = ax1.Position(1) + ax1.Position(3) - 0.2175;
cbh2.Position(1) = ax2.Position(1) + ax2.Position(3) - 0.2175;
cbh5.Position(1) = ax5.Position(1) + ax5.Position(3) - 0.2175;

% make colorbars level as subplots
cbh1.Position(2) = ax1.Position(2);
cbh2.Position(2) = ax2.Position(2);
cbh5.Position(2) = ax5.Position(2);

% DONE SETTING SUBPLOT POSITION
%--------------------------------------------------------

end   % end if pn == 33

% END MAKING y-avg. total tke PLOT
%--------------------------------------------------------------------------------

% (III_D) MAKE profile of hori. mean rho* total tke PLOT

if pn == 34 || pn == 20 
    
% find indices of x=-3 and x=3 (consistent with contours plot) in xh2, do the same for other cases
inxmin1 = find(abs(xh1-xmina) < 0.5*(xh1(2)-xh1(1)));   % index of xmin for dx=125
inxmin2 = find(abs(xh2-xmina) < 0.5*(xh2(2)-xh2(1)));   % index of xmin for dx=250 
inxmin5 = find(abs(xh5-xmina) < 0.5*(xh5(2)-xh5(1)));   % index of xmin for dx=500

inxmax1 = find(abs(xh1-xmaxa) < 0.5*(xh1(2)-xh1(1)));   % index of xmax for dx=125
inxmax2 = find(abs(xh2-xmaxa) < 0.5*(xh2(2)-xh2(1)));   % index of xmax for dx=250 
inxmax5 = find(abs(xh5-xmaxa) < 0.5*(xh5(2)-xh5(1)));   % index of xmax for dx=500

% get hori. mean (mean of y-avg.) total tke within xmin and xmax 
hmrhottke1 = mean(myrhottkem1r(inxmin1:inxmax1,:),1);   % hori. mean of total tke within region spec. by xmin and xmax, for dx=125
hmrhottke2 = mean(myrhottkem2r(inxmin2:inxmax2,:),1);   % hori. mean of total tke within region spec. by xmin and xmax, for dx=250
hmrhottke5 = mean(myrhottkem5r(inxmin5:inxmax5,:),1);   % hori. mean of total tke within region spec. by xmin and xmax, for dx=500

f = figure
% f = figure('units','normalized','outerposition',[0 0 1 1]);

% make profiles of hori. mean total tke for all cases
plot(hmrhottke1,z2d1(1,:),'-g','LineWidth',1.125,'color',[0 0.5 0]); hold on   % w-budget output
plot(hmrhottke2,z2d2(1,:),'-b','LineWidth',1.125); hold on   % w-budget output
plot(hmrhottke5,z2d5(1,:),'-r','LineWidth',1.125); hold on   % w-budget output

grid on;
box on;

lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m');
% lgd = legend('\(\Delta=62.5\) m','\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m','\(\Delta=1000\) m');
set(lgd,'Interpreter','latex','Location','northeast');

xlabel('\(\overline{\textrm{total TKE}}\) \big[kg m\textsuperscript{-1} s\textsuperscript{-2}\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

% format axes
ymin = 0;
ymax = 3;  
ytick = 1;

xmin = 0;
xmax = 1.5;
xtick = 0.5;

ax1 = gca;
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];
% ax1.FontSize = axfs;

pbr = pbra;
% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% set xlabel position (lower it)
ax1.XLabel.Position(2) = ax1.XLabel.Position(2) - 0.02;



end   % end if pn == 34
%---------------------------------------
% END MAKING profile of hori. mean total tke PLOT
%--------------------------------------------------------------------------------



% save plots
if pn == 11
filenameg = sprintf('myuwextke_u_%d%d',ti,tf);   % y-avg. u
end
if pn == 12
filenameg = sprintf('myuwextke_w_%d%d',ti,tf);   % y-avg. w
end
if pn == 2
filenameg = sprintf('myuwextke_ex_%d%d',ti,tf);   % y-avg. exn-star
end
if pn == 31
filenameg = sprintf('myuwextke_sgtke_%d%d',ti,tf);   % y-avg. subgrid tke
end
if pn == 32
filenameg = sprintf('myuwextke_retke_%d%d',ti,tf);   % y-avg. resolved tke
end
if pn == 33
filenameg = sprintf('myuwextke_ttke_%d%d',ti,tf);   % y-avg. total tke
end
if pn == 34
filenameg = sprintf('myuwextke_profile_ttke_%d%d',ti,tf);   % profile of hori. mean total tke
end

% saveas(f,fullfile(fnameg,filenameg),'epsc');



end   % end seq = 6:6
