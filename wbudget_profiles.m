close all
clear all

%----------------------------------------------------------------------------
% REMARKS:
% This code requires data generated from wbudget.m
% It produces profile of each term (hourly avg. of avg. along y 
% and within central 1 km in x)in w-budget model output for the cases of
% (dx=125,250, and 500 m)
%-----------------------------------------------------------------------------

tic 

% directory for dwdt_eqn.mat files
% fname = '/aos/home/mtang/Documents/matrices';
% fname = 'D:\Users\stang33\Desktop';
fname = 'C:\Users\SiuLung\Downloads';   % for dx = 250

% fielnames for wbudget.m 
filename1 = 'wbudget_prcl_dx125_dec03.mat';
filename2 = 'wbudget_prcl_dx250_dec03.mat';
filename5 = 'wbudget_prcl_dx500_dec03.mat';

% directory for output graphs
% fnameg = '/aos/home/mtang/Documents/graphs';
% fnameg = 'D:\Users\stang33\Desktop';
fnameg = 'C:\Users\SiuLung\Downloads';

% input range of z interested
zmax = 4;
zmin = 0;

% input plot number interested

% (1) = profile of each term in w-budget (1 panel)
% (2) = model rhs and dwdt (2 panels)

pn = 2;

%---------------------------------------------------------------------------------
% % get wbudget matrcies from wbudget.m
wbudget1 = load(fullfile(fname,filename1));  % terms array  for dx=125
wbudget2 = load(fullfile(fname,filename2));  % terms array  for dx=250
wbudget5 = load(fullfile(fname,filename5));  % terms array  for dx=500
% wbudget6 = load(fullfile(fname,filename6));  % terms array  for dx=62.5
% wbudget10 = load(fullfile(fname,filename10));  % terms array  for dx=1000

%----------------------------------------
for seq = 6:6
%    
ti = 0+(seq-1)*1; tf = 1+(seq-1)*1;
% 
% for P = 1:3   % subplot position
% 
fieldname = sprintf('f%d%d',ti,tf);

% get xh and z matrices for all plots 
zf = wbudget2.terms.(fieldname).zf;  % in km
xh1 = wbudget1.terms.(fieldname).xh;   % in km for dx=125
xh2 = wbudget2.terms.(fieldname).xh;   % in km for dx=250
xh5 = wbudget5.terms.(fieldname).xh;   % in km for dx=500

% get index of zmax on zf 
zfi = find(abs(zf-zmax)< 0.1/2);
zii = 1;

%-----------------------------------------------------
% (II) GET MATRICES FROM ARRAYS 

% get hadvm matrices from arrays (FROM W-BUDGET OUTPUT)
hadvm1 = wbudget1.terms.(fieldname).hadvm;   % time avg. of y-avg. hori. adv. for dx=125
hadvm2 = wbudget2.terms.(fieldname).hadvm;   % time avg. of y-avg. hori. adv. for dx=250
hadvm5 = wbudget5.terms.(fieldname).hadvm;   % time avg. of y-avg. hori. adv. for dx=500
% hadvm10 = wbudget10.terms.(fieldname).hadvm;   % time avg. of y-avg. hori. adv. for dx=1000
% hadvm6 = wbudget6.terms.(fieldname).hadvm;   % time avg. of y-avg. hori. adv. for dx=62.5

% reduce hadvm matrices into chosen region
hadvm1r = hadvm1(:,zii:zfi);   % for dx=125
hadvm2r = hadvm2(:,zii:zfi);   % for dx=250
hadvm5r = hadvm5(:,zii:zfi);   % for dx=500
% hadvm10r = hadvm10(:,zii:zfi);   % for dx=1000
% hadvm6r = hadvm6(:,zii:zfi);   % for dx=62.5

% get vadvm matrices from arrays (FROM W-BUDGET)
vadvm1 = wbudget1.terms.(fieldname).vadvm;   % time avg. of y-avg. vert. adv. for dx=125
vadvm2 = wbudget2.terms.(fieldname).vadvm;   % time avg. of y-avg. vert. adv. for dx=250
vadvm5 = wbudget5.terms.(fieldname).vadvm;   % time avg. of y-avg. vert. adv. for dx=500
% vadvm10 = wbudget10.terms.(fieldname).vadvm;   % time avg. of y-avg. vert. adv. for dx=1000
% vadvm6 = wbudget6.terms.(fieldname).vadvm;   % time avg. of y-avg. vert. adv. for dx=62.5

% reduce vadvm matrices into chosen region
vadvm1r = vadvm1(:,zii:zfi);   % for dx=125
vadvm2r = vadvm2(:,zii:zfi);   % for dx=250
vadvm5r = vadvm5(:,zii:zfi);   % for dx=500
% vadvm10r = vadvm10(:,zii:zfi);   % for dx=1000
% vadvm6r = vadvm6(:,zii:zfi);   % for dx=62.5

% get hidiffm matrices from arrays (FROM W-BUDGET OUTPUT)
hidiffm1 = wbudget1.terms.(fieldname).hidiffm;   % time avg. of y-avg. hori. diffusion for dx=125
hidiffm2 = wbudget2.terms.(fieldname).hidiffm;   % time avg. of y-avg. hori. diffusion for dx=250
hidiffm5 = wbudget5.terms.(fieldname).hidiffm;   % time avg. of y-avg. hori. diffusion for dx=500
% hidiffm10 = wbudget10.terms.(fieldname).hidiffm;   % time avg. of y-avg. hori. diffusion for dx=1000
% hidiffm6 = wbudget6.terms.(fieldname).hidiffm;   % time avg. of y-avg. hori. diffusion for dx=62.5

% reduce hidiffm matrices into chosen region
hidiffm1r = hidiffm1(:,zii:zfi);   % for dx=125
hidiffm2r = hidiffm2(:,zii:zfi);   % for dx=250
hidiffm5r = hidiffm5(:,zii:zfi);   % for dx=500
% hidiffm10r = hidiffm10(:,zii:zfi);   % for dx=1000
% hidiffm6r = hidiffm6(:,zii:zfi);   % for dx=62.5

% get vidiffm matrices from arrays (FROM W-BUDGET OUTPUT)
vidiffm1 = wbudget1.terms.(fieldname).vidiffm;   % time avg. of y-avg. vert. diffusion for dx=125
vidiffm2 = wbudget2.terms.(fieldname).vidiffm;   % time avg. of y-avg. vert. diffusion for dx=250
vidiffm5 = wbudget5.terms.(fieldname).vidiffm;   % time avg. of y-avg. vert. diffusion for dx=500
% vidiffm10 = wbudget10.terms.(fieldname).vidiffm;   % time avg. of y-avg. vert. diffusion for dx=1000
% vidiffm6 = wbudget6.terms.(fieldname).vidiffm;   % time avg. of y-avg. vert. diffusion for dx=62.5

% reduce vidiffm matrices into chosen region
vidiffm1r = vidiffm1(:,zii:zfi);   % for dx=125
vidiffm2r = vidiffm2(:,zii:zfi);   % for dx=250
vidiffm5r = vidiffm5(:,zii:zfi);   % for dx=500
% vidiffm10r = vidiffm10(:,zii:zfi);   % for dx=1000
% vidiffm6r = vidiffm6(:,zii:zfi);   % for dx=62.5


% get hturbm matrices from arrays (FROM W-BUDGET OUTPUT)
hturbm1 = wbudget1.terms.(fieldname).hturbm;   % time avg. of y-avg. hori. turb. for dx=125
hturbm2 = wbudget2.terms.(fieldname).hturbm;   % time avg. of y-avg. hori. turb. for dx=250
hturbm5 = wbudget5.terms.(fieldname).hturbm;   % time avg. of y-avg. hori. turb. for dx=500
% hturbm10 = wbudget10.terms.(fieldname).hturbm;   % time avg. of y-avg. hori. turb. for dx=1000
% hturbm6 = wbudget6.terms.(fieldname).hturbm;   % time avg. of y-avg. hori. turb. for dx=62.5

% reduce hturbm matrices into chosen region
hturbm1r = hturbm1(:,zii:zfi);   % for dx=125
hturbm2r = hturbm2(:,zii:zfi);   % for dx=250
hturbm5r = hturbm5(:,zii:zfi);   % for dx=500
% hturbm10r = hturbm10(:,zii:zfi);   % for dx=1000
% hturbm6r = hturbm6(:,zii:zfi);   % for dx=62.5

% get vturbm matrices from arrays (FROM W-BUDGET OUTPUT)
vturbm1 = wbudget1.terms.(fieldname).vturbm;   % time avg. of y-avg. vert. turb. for dx=125
vturbm2 = wbudget2.terms.(fieldname).vturbm;   % time avg. of y-avg. vert. turb. for dx=250
vturbm5 = wbudget5.terms.(fieldname).vturbm;   % time avg. of y-avg. vert. turb. for dx=500
% vturbm10 = wbudget10.terms.(fieldname).vturbm;   % time avg. of y-avg. vert. turb. for dx=1000
% vturbm6 = wbudget6.terms.(fieldname).vturbm;   % time avg. of y-avg. vert. turb. for dx=62.5

% reduce vturbm matrices into chosen region
vturbm1r = vturbm1(:,zii:zfi);   % for dx=125
vturbm2r = vturbm2(:,zii:zfi);   % for dx=250
vturbm5r = vturbm5(:,zii:zfi);   % for dx=500
% vturbm10r = vturbm10(:,zii:zfi);   % for dx=1000
% vturbm6r = vturbm6(:,zii:zfi);   % for dx=62.5

% get pgradm matrices from arrays (FROM W-BUDGET OUTPUT)
pgradm1 = wbudget1.terms.(fieldname).pgradm;   % time avg. of y-avg. pgf for dx=125
pgradm2 = wbudget2.terms.(fieldname).pgradm;   % time avg. of y-avg. pgf for dx=250
pgradm5 = wbudget5.terms.(fieldname).pgradm;   % time avg. of y-avg. pgf for dx=500
% pgradm10 = wbudget10.terms.(fieldname).pgradm;   % time avg. of y-avg. pgf for dx=1000
% pgradm6 = wbudget6.terms.(fieldname).pgradm;   % time avg. of y-avg. pgf for dx=62.5

% reduce pgradm matrices into chosen region
pgradm1r = pgradm1(:,zii:zfi);   % for dx=125
pgradm2r = pgradm2(:,zii:zfi);   % for dx=250
pgradm5r = pgradm5(:,zii:zfi);   % for dx=500
% pgradm10r = pgradm10(:,zii:zfi);   % for dx=1000
% pgradm6r = pgradm6(:,zii:zfi);   % for dx=62.5

% get rdampm matrices from arrays (FROM W-BUDGET OUTPUT)
rdampm1 = wbudget1.terms.(fieldname).rdampm;   % time avg. of y-avg. rdampm for dx=125
rdampm2 = wbudget1.terms.(fieldname).rdampm;   % time avg. of y-avg. rdampm for dx=250
rdampm5 = wbudget1.terms.(fieldname).rdampm;   % time avg. of y-avg. rdampm for dx=500

% reduce rdampm matrices into chosen region
rdampm1r = rdampm1(:,zii:zfi);   % for dx=125
rdampm2r = rdampm2(:,zii:zfi);   % for dx=250
rdampm5r = rdampm5(:,zii:zfi);   % for dx=500

% get buoym matrices from arrays (FROM W-BUDGET)
buoym1w = wbudget1.terms.(fieldname).buoym;   % time avg. of y-avg. buoym for dx=125
buoym2w = wbudget2.terms.(fieldname).buoym;   % time avg. of y-avg. buoym for dx=250
buoym5w = wbudget5.terms.(fieldname).buoym;   % time avg. of y-avg. buoym for dx=500
% buoym10w = wbudget10.terms.(fieldname).buoym;   % time avg. of y-avg. buoym for dx=1000
% buoym6w = wbudget6.terms.(fieldname).buoym;   % time avg. of y-avg. buoym for dx=62.5

% reduce buoymw matrices into chosen region
buoym1rw = buoym1w(:,zii:zfi);   % for dx=125
buoym2rw = buoym2w(:,zii:zfi);   % for dx=250
buoym5rw = buoym5w(:,zii:zfi);   % for dx=500
% buoym10rw = buoym10w(:,zii:zfi);   % for dx=1000
% buoym6rw = buoym6w(:,zii:zfi);   % for dx=62.5

% get myrhsm matrices from wbudget arrays (FROM W-BUDGET OUTPUT)

myrhsm1w = wbudget1.terms.(fieldname).myrhsm;   % time avg. of y-avg. r.h.s. terms for dx=125
myrhsm2w = wbudget2.terms.(fieldname).myrhsm;   % time avg. of y-avg. r.h.s. terms for dx=250
myrhsm5w = wbudget5.terms.(fieldname).myrhsm;   % time avg. of y-avg. r.h.s. terms for dx=500
% myrhsm10w = wbudget10.terms.(fieldname).myrhsm;   % time avg. of y-avg. r.h.s. terms for dx=1000
% myrhsm6w = wbudget6.terms.(fieldname).myrhsm;   % time avg. of y-avg. r.h.s. terms for dx=62.5

% reduce myrhsm matrices into chosen region
myrhsm1r = myrhsm1w(:,zii:zfi);   % for dx=125
myrhsm2r = myrhsm2w(:,zii:zfi);   % for dx=250
myrhsm5r = myrhsm5w(:,zii:zfi);   % for dx=500
% myrhsm10r = myrhsm10w(:,zii:zfi);   % for dx=1000
% myrhsm6r = myrhsm6w(:,zii:zfi);   % for dx=62.5

% get ddtwm matrices from arrays (FROM W-BUDGET OUTPUT)
ddtwm1w = wbudget1.terms.(fieldname).ddtwm;   % time avg. of y-avg. dwdt for dx=125
ddtwm2w = wbudget2.terms.(fieldname).ddtwm;   % time avg. of y-avg. dwdt for dx=250
ddtwm5w = wbudget5.terms.(fieldname).ddtwm;   % time avg. of y-avg. dwdt for dx=500
% ddtwm10w = wbudget10.terms.(fieldname).ddtwm;   % time avg. of y-avg. dwdt for dx=1000
% ddtwm6w = wbudget6.terms.(fieldname).ddtwm;   % time avg. of y-avg. dwdt for dx=62.5

% reduce ddtwm matrices into chosen region
ddtwm1rw = ddtwm1w(:,zii:zfi);   % for dx=125
ddtwm2rw = ddtwm2w(:,zii:zfi);   % for dx=250
ddtwm5rw = ddtwm5w(:,zii:zfi);   % for dx=500
% ddtwm10rw = ddtwm10w(:,zii:zfi);   % for dx=1000
% ddtwm6rw = ddtwm6w(:,zii:zfi);   % for dx=62.5

%-----------------------------------------------------

% (III) GET MEAN OVER CENTRAL 1 KM OF X-DOMAIN

% NOTE: both dwdt_eqn_2.m and wbudget.m give variables in half-level xh

% for dx=125
dxh1 = xh1(2:end) - xh1(1:end-1);   % dxh for dx=125 case
xmin1 = find(abs(xh1-(-0.5)) < dxh1(1)/2);   % index of -0.5 on xh1
xmax1 = find(abs(xh1-(0.5)) < dxh1(1)/2);   % index of +0.5 on xh1

% for dx=250
dxh2 = xh2(2:end) - xh2(1:end-1);   % dxh for dx=250 case
xmin2 = find(abs(xh2-(-0.5)) < dxh2(1)/2);   % index of -0.5 on xh2
xmax2 = find(abs(xh2-(0.5)) < dxh2(1)/2);   % index of +0.5 on xh2

% for dx=500
dxh5 = xh5(2:end) - xh5(1:end-1);   % dxh for dx=500 case
xmin5 = max(find(abs(xh5-(-0.5)) < dxh5(1)/1.5));   % index of -0.5 on xh5
xmax5 = min(find(abs(xh5-(0.5)) < dxh5(1)/1.5));   % index of +0.5 on xh5
%-------------------------------------------

% get mean y-avg. hori. adv. (FROM W-BUDGET)
% hadvm6rm = mean(hadvm6r(xmin6:xmax6,:),1);   % for dx=62.5
hadvm1rm = mean(hadvm1r(xmin1:xmax1,:),1);   % for dx=125
hadvm2rm = mean(hadvm2r(xmin2:xmax2,:),1);   % for dx=250
hadvm5rm = mean(hadvm5r(xmin5:xmax5,:),1);   % for dx=500
% hadvm10rm = mean(hadvm10r(xmin10:xmax10,:),1);   % for dx=1000

% get mean y-avg. vert. adv. (FROM W-BUDGET)
% vadvm6rm = mean(vadvm6r(xmin6:xmax6,:),1);   % for dx=62.5
vadvm1rm = mean(vadvm1r(xmin1:xmax1,:),1);   % for dx=125
vadvm2rm = mean(vadvm2r(xmin2:xmax2,:),1);   % for dx=250
vadvm5rm = mean(vadvm5r(xmin5:xmax5,:),1);   % for dx=500
% vadvm10rm = mean(vadvm10r(xmin10:xmax10,:),1);   % for dx=1000

% get mean y-avg. hori. diff. (FROM W-BUDGET)
% hidiffm6rm = mean(hidiffm6r(xmin6:xmax6,:),1);   % for dx=62.5
hidiffm1rm = mean(hidiffm1r(xmin1:xmax1,:),1);   % for dx=125
hidiffm2rm = mean(hidiffm2r(xmin2:xmax2,:),1);   % for dx=250
hidiffm5rm = mean(hidiffm5r(xmin5:xmax5,:),1);   % for dx=500
% hidiffm10rm = mean(hidiffm10r(xmin10:xmax10,:),1);   % for dx=1000

% get mean y-avg. vert. diff. (FROM W-BUDGET)
% vidiffm6rm = mean(vidiffm6r(xmin6:xmax6,:),1);   % for dx=62.5
vidiffm1rm = mean(vidiffm1r(xmin1:xmax1,:),1);   % for dx=125
vidiffm2rm = mean(vidiffm2r(xmin2:xmax2,:),1);   % for dx=250
vidiffm5rm = mean(vidiffm5r(xmin5:xmax5,:),1);   % for dx=500
% vidiffm10rm = mean(vidiffm10r(xmin10:xmax10,:),1);   % for dx=1000

% get mean y-avg. hori. turb. (FROM W-BUDGET)
% hturbm6rm = mean(hturbm6r(xmin6:xmax6,:),1);   % for dx=62.5
hturbm1rm = mean(hturbm1r(xmin1:xmax1,:),1);   % for dx=125
hturbm2rm = mean(hturbm2r(xmin2:xmax2,:),1);   % for dx=250
hturbm5rm = mean(hturbm5r(xmin5:xmax5,:),1);   % for dx=500
% hturbm10rm = mean(hturbm10r(xmin10:xmax10,:),1);   % for dx=1000

% get mean y-avg. vert. turb. (FROM W-BUDGET)
% vturbm6rm = mean(vturbm6r(xmin6:xmax6,:),1);   % for dx=62.5
vturbm1rm = mean(vturbm1r(xmin1:xmax1,:),1);   % for dx=125
vturbm2rm = mean(vturbm2r(xmin2:xmax2,:),1);   % for dx=250
vturbm5rm = mean(vturbm5r(xmin5:xmax5,:),1);   % for dx=500
% vturbm10rm = mean(vturbm10r(xmin10:xmax10,:),1);   % for dx=1000

% get mean y-avg. pgf (FROM W-BUDGET)
% pgradm6rm = mean(pgradm6r(xmin6:xmax6,:),1);   % for dx=62.5
pgradm1rm = mean(pgradm1r(xmin1:xmax1,:),1);   % for dx=125
pgradm2rm = mean(pgradm2r(xmin2:xmax2,:),1);   % for dx=250
pgradm5rm = mean(pgradm5r(xmin5:xmax5,:),1);   % for dx=500
% pgradm10rm = mean(pgradm10r(xmin10:xmax10,:),1);   % for dx=1000

% get mean y-avg. rdamp (FROM W-BUDGET)
rdampm1rm = mean(rdampm1r(xmin1:xmax1,:),1);   % for dx=125
rdampm2rm = mean(rdampm2r(xmin2:xmax2,:),1);   % for dx=250
rdampm5rm = mean(rdampm5r(xmin5:xmax5,:),1);   % for dx=500

% get mean y-avg. buoy. (FROM W-BUDGET)
% buoym6rmw = mean(buoym6rw(xmin6:xmax6,:),1);   % for dx=62.5
buoym1rmw = mean(buoym1rw(xmin1:xmax1,:),1);   % for dx=125
buoym2rmw = mean(buoym2rw(xmin2:xmax2,:),1);   % for dx=250
buoym5rmw = mean(buoym5rw(xmin5:xmax5,:),1);   % for dx=500
% buoym10rmw = mean(buoym10rw(xmin10:xmax10,:),1);   % for dx=1000

%-----------------------------------------
% extrapolate buoy. profile to the surface
bslope = 0.1/(buoym2rmw(3)-buoym2rmw(2));
bsfc = buoym2rmw(2) -0.1/bslope;
buoym2rmw(1) = bsfc;
%-----------------------------------------

% get mean y-avg. r.h.s. (FROM W-BUDGET OUTPUT)
% myrhsm6rm = mean(myrhsm6r(xmin6:xmax6,:),1);   % for dx=62.5
myrhsm1rm = mean(myrhsm1r(xmin1:xmax1,:),1);   % for dx=125
myrhsm2rm = mean(myrhsm2r(xmin2:xmax2,:),1);   % for dx=250
myrhsm5rm = mean(myrhsm5r(xmin5:xmax5,:),1);   % for dx=500
% myrhsm10rm = mean(myrhsm10r(xmin10:xmax10,:),1);   % for dx=1000





% get mean y-avg. dwdt (FROM W-BUDGET OUTPUT)
% ddtw6rwm = mean(ddtwm6rw(xmin6:xmax6,:),1);   % for dx=62.5
ddtw1rwm = mean(ddtwm1rw(xmin1:xmax1,:),1);   % for dx=125
ddtw2rwm = mean(ddtwm2rw(xmin2:xmax2,:),1);   % for dx=250
ddtw5rwm = mean(ddtwm5rw(xmin5:xmax5,:),1);   % for dx=500
% ddtw10rwm = mean(ddtwm10rw(xmin10:xmax10,:),1);   % for dx=1000
%-----------------------------------------------
% SET PLOT ATTRIBUTES FONT SIZE

tfs = 16;   % title font size
tefs = 18;   % text font size (plot number)
axfs = 14;   % axes font size

% if pn == 1
% tfs = 13;   % title font size
% tefs = 13;   % text font size (plot number)
% axfs = 11;   % axes font size
% end

pbra = 2;   % apsect ratio for plots, y:x ratio


%-----------------------------------------------------

% (IV) MALE PROFILE FOR EACH TERM IN W-BUDGET

if pn == 1
    
f = figure('units','normalized','outerposition',[0 0 1 1])
% f = figure

%----------------------------------------------------------------------
% (V_A) MAKE PLOT OF PROFILE OF R.H.S. OF W-BUDGET

% % the 1st entry in myrhsm matrices is an outlier, set it equal to the 2nd
% % entry
% myrhsm1rmnew = [myrhsm1rm(2) myrhsm1rm(1,2:end)];   % for dx=125
% myrhsm2rmnew = [myrhsm2rm(1,2) myrhsm2rm(1,2:end)];   % for dx=250
% myrhsm5rmnew = [myrhsm5rm(2) myrhsm5rm(1,2:end)];   % for dx=500

% plot(hadvm2rm,zf(zii:zfi),'-.b','LineWidth',1.125,'color',[0, 0.4470, 0.7410]);   hold on   % w-budget output, light blue
% plot(vadvm2rm,zf(zii:zfi),'--b','LineWidth',1.125,'color',[0, 0.4470, 0.7410]);   hold on   % w-budget output
plot(hadvm2rm,zf(zii:zfi),'-.b','LineWidth',1.5,'color',[0, 0.4470, 0.7410]);   hold on   % w-budget output, light blue
plot(vadvm2rm,zf(zii:zfi),'--b','LineWidth',1.5,'color',[0, 0.4470, 0.7410]);   hold on   % w-budget output

% plot(hidiffm2rm,zf(zii:zfi),'-.b','LineWidth',1.125,'color',[0.6350, 0.0780, 0.1840]);   hold on   % w-budget output, maroon
% plot(vidiffm2rm,zf(zii:zfi),'--b','LineWidth',1.125,'color',[0.6350, 0.0780, 0.1840]);   hold on   % w-budget output

plot(hidiffm2rm,zf(zii:zfi),'-.b','LineWidth',1.5,'color',[0.6350, 0.0780, 0.1840]);   hold on   % w-budget output, maroon
plot(vidiffm2rm,zf(zii:zfi),'--b','LineWidth',1.5,'color',[0.6350, 0.0780, 0.1840]);   hold on   % w-budget output

% plot(hturbm2rm,zf(zii:zfi),'-.b','LineWidth',1.125,'color',[0.25, 0.25, 0.25]);   hold on   % w-budget output, grey
% plot(vturbm2rm,zf(zii:zfi),'--b','LineWidth',1.125,'color',[0.25, 0.25, 0.25]);   hold on   % w-budget output

plot(hturbm2rm,zf(zii:zfi),'-.b','LineWidth',1.5,'color',[0.25, 0.25, 0.25]);   hold on   % w-budget output, grey
plot(vturbm2rm,zf(zii:zfi),'--b','LineWidth',1.5,'color',[0.25, 0.25, 0.25]);   hold on   % w-budget output


% plot(pgradm1rm,zf(zii:zfi),'--g','LineWidth',1.125);   hold on   % w-budget output, blue
plot(pgradm2rm,zf(zii:zfi),'--b','LineWidth',1.5,'color',[0 0 0.5]);   hold on   % w-budget output, blue
% plot(pgradm5rm,zf(zii:zfi),'--r','LineWidth',1.125);   hold on   % w-budget output, blue


% plot(buoym1rmw,zf(zii:zfi),'-b','LineWidth',1.125,'color',[0 0.5 0]);   hold on  % w-budget output
plot(buoym2rmw,zf(zii:zfi),'-b','LineWidth',1.5,'color',[0 0 0.5]);   hold on  % w-budget output
% plot(buoym5rmw,zf(zii:zfi),'-r','LineWidth',1.125);   hold off  % w-budget output

plot(rdampm2rm,zf(zii:zfi),'-.','LineWidth',1.5,'color',[0 0 0.5]);   hold on  % w-budget output

% myrhsm2rmnew = [myrhsm2rm(1,2) myrhsm2rm(1,2:end)];   % for dx=250
% plot(myrhsm2rmnew,zf(zii:zfi),'--g','LineWidth',1.125,'color',[0 0 0]); hold on   % w-budget output

grid on;
box on;

lgd = legend('hori. adv.','vert. adv.','hori. diff.','vert. diff','hori. turb.','vert. turb.','pgf','buoy.','damping');
% % lgd = legend('\(\Delta=62.5\) m','\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m','\(\Delta=1000\) m');
% set(lgd,'Interpreter','latex','Location','northeast');
% set(lgd,'Units','normalized','Interpreter','latex','Location','northeastoutside');
set(lgd,'Units','normalized','Interpreter','latex');
lgd.Position = [0.635 0.684 0.07 0.2];
% % 
xlabel('\big[ms\textsuperscript{-2}\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['Model R.H.S., ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% % title with no time info.
% title(['R.H.S.'],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 3;   % originally 3 
ytick = 1;   % originally 1

xmin = -0.08;
xmax = 0.08;
xtick = 0.04;

if seq == 5
xmin = -2e-4;
xmax = 2e-4;
xtick = 2e-4/2; 
end

% get current axes and store it to variable named ax1
ax1 = gca;
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];
ax1.FontSize = axfs;

pbr = pbra;
% % Make y-axis:x-axis = 3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(0.95,1.025,'(a)','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% text(0.63,0.75,'imoist\(=0\)','Units','normalized','Interpreter', 'latex','FontSize',tefs-3,'FontWeight','bold');
% text(0.58,0.72,'half-domain','Units','normalized','Interpreter', 'latex','FontSize',tefs-3,'FontWeight','bold');

% text(0.17,0.69,'avg. over \(-0.5\leq x\leq0.5\)','Units','normalized','Interpreter', 'latex','FontSize',tefs-3,'FontWeight','bold');

end   % end if pn == 1
% MAKING PLOT FOR EACH TERM IN W-BUDGET
%--------------------------------------------------------------


% (V) MALE PROFILE FOR R.H.S. and dw/dt FOR W-BUDGET

if pn == 2
    
f = figure('units','normalized','outerposition',[0 0 1 1])
% f = figure

%----------------------------------------------------------------------
% (V_A) MAKE PLOT OF PROFILE OF R.H.S. OF W-BUDGET

subplot(1,2,1)

% the 1st entry in myrhsm matrices is an outlier, set it equal to the 2nd
% entry
myrhsm1rmnew = myrhsm1rm;   % for dx=125
myrhsm2rmnew = myrhsm2rm;   % for dx=250
myrhsm5rmnew = myrhsm5rm;   % for dx=500

% myrhsm1rmnew = [myrhsm1rm(2) myrhsm1rm(1,2:end)];   % for dx=125
% myrhsm2rmnew = [myrhsm2rm(1,2) myrhsm2rm(1,2:end)];   % for dx=250
% myrhsm5rmnew = [myrhsm5rm(2) myrhsm5rm(1,2:end)];   % for dx=500

% plot(myrhsm6rm,zf(zii:zfi),'--m','LineWidth',1.125); hold on   % w-budget output
plot(myrhsm1rmnew,zf(zii:zfi),'--g','LineWidth',1.125,'color',[0 0.5 0]); hold on   % w-budget output
plot(myrhsm2rmnew,zf(zii:zfi),'--b','LineWidth',1.125); hold on   % w-budget output
plot(myrhsm5rmnew,zf(zii:zfi),'--r','LineWidth',1.125); hold on   % w-budget output
% plot(myrhsm10rm,zf(zii:zfi),'--','LineWidth',1.125,'color',[1 0.5 0]); hold off   % w-budget output

grid on;
box on;

%lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m');
% % lgd = legend('\(\Delta=62.5\) m','\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m','\(\Delta=1000\) m');
%set(lgd,'Interpreter','latex','Location','northeast');

xlabel('\big[ms\textsuperscript{-2}\big]','Interpreter','latex','FontSize',12);
ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['Model R.H.S., ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% title with no time info.
title(['R.H.S.'],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 4;  
ytick = 1;

xmin = -1e-4;
xmax = 1e-4;
xtick = 1e-4/2;

ax1 = gca;
ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

ax1.YLim = [ymin,ymax];
ax1.YTick = [ymin:ytick:ymax];
ax1.FontSize = axfs;

pbr = pbra;
% % Make y-axis:x-axis = 3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(0.95,1.03,'(a)','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% text(0.63,0.75,'imoist\(=0\)','Units','normalized','Interpreter', 'latex','FontSize',tefs-3,'FontWeight','bold');
% text(0.58,0.72,'half-domain','Units','normalized','Interpreter', 'latex','FontSize',tefs-3,'FontWeight','bold');

% text(0.17,0.69,'avg. over \(-0.5\leq x\leq0.5\)','Units','normalized','Interpreter', 'latex','FontSize',tefs-3,'FontWeight','bold');

% MAKING PLOT FOR R.H.S. OF W-BUDGET
%--------------------------------------------------------------

% (V_B) MAKE PLOT OF PROFILE FOR dw/dt of W-BUDGET

subplot(1,2,2)
% plot(ddtw6rwm,zf(zii:zfi),'-m','LineWidth',1.125); hold on   % w-budget output
plot(ddtw1rwm,zf(zii:zfi),'-g','LineWidth',1.125,'color',[0 0.5 0]); hold on   % w-budget output
plot(ddtw2rwm,zf(zii:zfi),'-b','LineWidth',1.125); hold on   % w-budget output
plot(ddtw5rwm,zf(zii:zfi),'-r','LineWidth',1.125); hold on   % w-budget output
% plot(ddtw10rwm,zf(zii:zfi),'-','LineWidth',1.125,'color',[1 0.5 0]); hold off   % w-budget output

grid on;
box on;

lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m');
% lgd = legend('\(\Delta=62.5\) m','\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m','\(\Delta=1000\) m');
set(lgd,'Interpreter','latex','Location','northeast');

xlabel('\big[ms\textsuperscript{-2}\big]','Interpreter','latex','FontSize',12);
% ylabel('\textit{z}-position \big[km\big]','Interpreter','latex','FontSize',12);

timehi = sprintf('%02d',ti);
timemini = sprintf('%02d',0);
timehf = sprintf('%02d',tf);
timeminf = sprintf('%02d',0);

% title(['\(d\overline{w}/dt\), ',num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf)...
%     ':' num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% title with no time info.
title(['\(d\overline{w}/dt\)'],'Interpreter','latex','FontSize',tfs);

% format axes
ymin = 0;
ymax = 4;  
ytick = 1;

xmin = -1e-4;
xmax = 1e-4;
xtick = 1e-4/2;

if seq == 5
xmin = -2e-4;
xmax = 2e-4;
xtick = 2e-4/2; 
end

ax2 = gca;
ax2.XLim = [xmin,xmax];
ax2.XTick = [xmin:xtick:xmax];

ax2.YLim = [ymin,ymax];
ax2.YTick = [ymin:ytick:ymax];
ax2.FontSize = axfs;

% set ax2 left position closer to ax1
ax2.Position(1) = ax1.Position(1) + 0.72*ax1.Position(3);

pbr = pbra;
% % Make y-axis:x-axis = 3
pbaspect([(xmax-xmin),(pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

text(0.95,1.03,'(b)','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% text(0.63,0.75,'imoist\(=0\)','Units','normalized','Interpreter', 'latex','FontSize',tefs-3,'FontWeight','bold');
% text(0.58,0.72,'half-domain','Units','normalized','Interpreter', 'latex','FontSize',tefs-3,'FontWeight','bold');

% text(0.17,0.69,'avg. over \(-0.5\leq x\leq0.5\)','Units','normalized','Interpreter', 'latex','FontSize',tefs-3,'FontWeight','bold');

% DONE MAKING PROFILE OF dw/dt OF W-BUDGET
%--------------------------------------------------------------

% DONE MAKING PLOT FOR R.H.S. AND dw/dt OF W-BUDGET
%--------------------------------------------------------------

end   % end if pn == 2
%--------------------------------------------------------------

% save plots
if pn == 1
filenameg = sprintf('wbudget_profile_inditerm_%d%d',ti,tf);   % comparison of individual term in dwdt eqb.
end
if pn == 2
filenameg = sprintf('wbudget_profile_rhsdwdt_%d%d',ti,tf);   % profiles of r.h.s. and dw/dt 
end

% saveas(f,fullfile(fnameg,filenameg),'epsc');




end   % end for seq = 6:6