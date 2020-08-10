close all
clear all

tic 

% directory for pdata.mat files
% fname = '/aos/home/mtang/Documents/matrices';
fname = 'C:\Users\SiuLung\Downloads\data_paper';   % for dx = 250
% filename1 = 'prcl_wlcl_dx125_aug02.mat';
% filename2 = 'prcl_wlcl_dx250_aug02.mat';
% filename5 = 'prcl_wlcl_dx500_aug02.mat';

filename1 = 'prcl_wlcl_dx125_aug02_1.mat';
filename2 = 'prcl_wlcl_dx250_aug02_1.mat';
filename5 = 'prcl_wlcl_dx500_aug02_1.mat';

% directory for output graphs
% fnameg = 'aos/home/mtang/Documents/graphs';
fnameg = 'C:\Users\SiuLung\Downloads';

%-------------------------------------------------------------

% % get pdata matrcies
pdata1 = load(fullfile(fname,filename1));  % pdata array  for dx=125
pdata2 = load(fullfile(fname,filename2));  % pdata array  for dx=250
pdata5 = load(fullfile(fname,filename5));  % pdata array  for dx=500

% % get prclzlcl matrix
prclzlcl1 = pdata1.pdata.prclzlcl;  % in m for dx=125 in m
prclzlcl2 = pdata2.pdata.prclzlcl;   % in m for dx=250 in m
prclzlcl5 = pdata5.pdata.prclzlcl;   % in m for dx=500 in m

% get total number of parcels for each case
numpr1 = size(prclzlcl1,1);   % number of parcels for dx=125
numpr2 = size(prclzlcl2,1);   % number of parcels for dx=250
numpr5 = size(prclzlcl5,1);   % number of parcels for dx=500

% get max. height reached by each parcel
maxh1 = max(prclzlcl1,[],2);   % max. height reached by each parcel for dx=125 in m
maxh2 = max(prclzlcl2,[],2);   % max. height reached by each parcel for dx=250 in m
maxh5 = max(prclzlcl5,[],2);   % max. height reached by each parcel for dx=500 in m

% get max. height in km
maxh1k = maxh1./1000;   % for dx=125
maxh2k = maxh2./1000;   % for dx=250
maxh5k = maxh5./1000;   % for dx=500

lper = 11;

% initialize percent matrices 
percent1 = zeros(1,lper);   % for dx=125
percent2 = zeros(1,lper);   % for dx=250
percent5 = zeros(1,lper);   % for dx=500

% for dx=125
for i = 1:lper
    if find(maxh1k >= 1+(i-1)*0.5) ~= 0
        num1 = length(find(maxh1k >= 1+(i-1)*0.5));
        percent1(i) = num1/numpr1;
    else 
        percent1(i) = 0;
    end
end

% for dx=250
for i = 1:lper
    if find(maxh2k >= 1+(i-1)*0.5) ~= 0
        num2 = length(find(maxh2k >= 1+(i-1)*0.5));
        percent2(i) = num2/numpr2;
    else 
        percent2(i) = 0;
    end
end

% for dx=500
for i = 1:lper
    if find(maxh5k >= 1+(i-1)*0.5) ~= 0
        num5 = length(find(maxh5k >= 1+(i-1)*0.5));
        percent5(i) = num5/numpr5;
    else 
        percent5(i) = 0;
    end
end

%--------------------------------------------------------------------------
% MAKE BAR PLOT

% centers of bins
xbins = [1:0.5:6]; 
% xbins = [0:0.5:6];

% concatenate percent matrices together
percentall = [percent1;percent2;percent5];

% take transpose such that each column represents pdf of each grid spacing
pdf = percentall';

% open figure
f = figure('units','normalized','outerposition',[0 0 1 1])
% figure

sp = subplot(2,1,2);
% b = bar(xbins,100.*pdf);
b = barh(xbins,100.*pdf,1.2,'hist');

grid on;
box on;

% set color for each grid spacing
b(1).FaceColor = [0 0.5 0];
b(2).FaceColor = 'b';
b(3).FaceColor = 'r';

% set legend
lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m');
set(lgd,'Interpreter','latex','Location','northeast');

% set x- and y-labels
xlabel('\%','Interpreter','latex','FontSize',12);
ylabel('Height \big[km\big]','Interpreter','latex','FontSize',12);

% format axis
ax1 = gca;
% ax.XTick = xbins;
% ax.YTick = xbins;

ymin = 0.5;
ymax = 6.5;  

xmin = 0;
xmax = 100;
xtick = 10;


ax1.YLim = [ymin,ymax];
ax1.YTick = xbins;

ax1.XLim = [xmin,xmax];
ax1.XTick = [xmin:xtick:xmax];

pbr = 2.5;
% % Make -axis:x-axis = 2
pbaspect([(xmax-xmin),(1/pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(0.975,1.075,sprintf('(a)'),'Units','normalized','Interpreter', 'latex','FontSize',15,'FontWeight','bold');
text(0.975,1.075,sprintf('(b)'),'Units','normalized','Interpreter', 'latex','FontSize',15,'FontWeight','bold');

%--------------------------------------------------------------------------
% MAKE WLCL PLOT

% % get wlcl matrix
wlcl1 = pdata1.pdata.wlclr;  % in m/s for dx=125 in m
wlcl2 = pdata2.pdata.wlclr;   % in m/s for dx=250 in m
wlcl5 = pdata5.pdata.wlclr;   % in m/s for dx=500 in m

% % get max. wlcl of each parcel
maxw1 = max(wlcl1,[],2);   % max. w by each parcel for dx=125 in m
maxw2 = max(wlcl2,[],2);   % max. w by each parcel for dx=250 in m
maxw5 = max(wlcl5,[],2);   % max. w by each parcel for dx=500 in m

% get length of each matrix
lwlcl1 = length(maxw1);  % for dx=125
lwlcl2 = length(maxw2);  % for dx=250
lwlcl5 = length(maxw5);  % for dx=500

% set edges of bins for bar plot
edges = 0:0.5:5.5;

% get counts of wlcl in each bin (no. of wlcl falls into certain bin)
counts1 = histcounts(wlcl1,edges);   % for dx=125
counts2 = histcounts(wlcl2,edges);   % for dx=250
counts5 = histcounts(wlcl5,edges);   % for dx=500

% get prob. of finding wlcl in certain bin 
prob1 = counts1./(sum(counts1));   % for dx=125
prob2 = counts2./(sum(counts2));   % for dx=250
prob5 = counts5./(sum(counts5));   % for dx=500

% centers of bins
wlclxbins = [0.25:0.5:5.25];
% xbins = [0:0.5:6];

% concatenate prob. matrices together
proball = [prob1;prob2;prob5];

% take transpose such that each column represents prob. of each grid spacing
% proballt = proball';
proballt = proball';

% get probability density function from proballt
pdfdata = proballt.*(1/(wlclxbins(2)-wlclxbins(1)))';

% open figure
% f2 = figure('units','normalized','outerposition',[0 0 1 1])
% figure

sp = subplot(2,1,1);
% b = bar(xbins,100.*pdf);
b = bar(wlclxbins,proballt,0.8,'hist'); 
% b = bar(wlclxbins,pdfdata,0.8,'hist'); hold on
% plot(wlclxbins,pdfdata(:,1),'Color',[0 0.5 0]); hold on
% plot(wlclxbins,pdfdata(:,2),'b'); hold on
% plot(wlclxbins,pdfdata(:,3),'r'); hold off
% y = poisspdf(wlclxbins, 0.09);

grid on;
box on;

% set color for each grid spacing
b(1).FaceColor = [0 0.5 0];
b(2).FaceColor = 'b';
b(3).FaceColor = 'r';

% set legend
lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m');
set(lgd,'Interpreter','latex','Location','northeast');

% set x- and y-labels
% xlabel('w\(_{LCL}\) \big[ms\textsuperscript{-1}\big]','Interpreter','latex','FontSize',12);
xlabel('\(w_{\mathrm{LCL}}\) \big[ms\textsuperscript{-1}\big]','Interpreter','latex','FontSize',12);
% ylabel('\%','Interpreter','latex','FontSize',12);
ylabel('\(N_{w}/N_{\mathrm{LCL}}\)','Interpreter','latex','FontSize',12);

% format axis
ax2 = gca;
% ax.XTick = xbins;
% ax.XTick = wlclxbins;

ymin = 0;
ymax = 0.5;  
ytick = 0.1;

xmin = 0;
xmax = 5.5;

ax2.YLim = [ymin,ymax];
ax2.YTick = [ymin:ytick:ymax];

ax2.XLim = [xmin,xmax];
ax2.XTick = edges;

pbr = 2.5;
% % Make -axis:x-axis = 2
pbaspect([(xmax-xmin),(1/pbr)*(xmax-xmin),1]); % multiple y-axis by the factor

% text(0.975,1.075,sprintf('(b)'),'Units','normalized','Interpreter', 'latex','FontSize',15,'FontWeight','bold');
text(0.975,1.075,sprintf('(a)'),'Units','normalized','Interpreter', 'latex','FontSize',15,'FontWeight','bold');

%--------------------------------------------------------------------------

% MAKE THE TWO SUBPLOTS CLOSER TO EACH OTHER
ax1.Position(2) = ax1.Position(2) + 0.05;


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% SAVE FIGURE

filenameg = sprintf('prcl_track_wlcl');
% saveas(f,fullfile(fnameg,filenameg),'epsc');




