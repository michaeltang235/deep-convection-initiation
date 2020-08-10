close all 
clear all

%----------------------------------------------------------------------------
% REMARKS:
% this code makes plots of y-avg. of time avg. sfc. precipitation as a
% function of x
%----------------------------------------------------------------------------
tic 

% directory for dwdt_eqn.mat files
% fname = '/aos/home/mtang/Documents/matrices';
fname = 'C:\Users\SiuLung\Downloads';   % for dx = 250

% enter filenames
filenamect1 = 'cth_dx125_may08_1.txt';
filenamel1 = 'lnb_dx125_may08.txt';

filenamect2 = 'cth_dx250_mar07_1.txt';
filenamel2 = 'lnb_dx250_mar07.txt';

filenamect5 = 'cth_dx500_feb27_1.txt';
filenamel5 = 'lnb_dx500_feb27.txt';

filenamect10 = 'cth_dx1000_nov08_1.txt';

% get cloud-top height matrices
cth125 = load(fullfile(fname,filenamect1));
lnb125 = load(fullfile(fname,filenamel1));

cth250 = load(fullfile(fname,filenamect2));
lnb250 = load(fullfile(fname,filenamel2));

cth500 = load(fullfile(fname,filenamect5));
lnb500 = load(fullfile(fname,filenamel5));

cth1000 = load(fullfile(fname,filenamect10));

%----------------------------------------------------------------------------
% MAKE PLOT of CLOUD-TOP HEIGHT

ts1 = ((1:1:121)-1).*360;  % make time series for 121 input files
ts2 = ((1:1:481)-1).*90;  % make time series for 481 input files

f = figure

p1 = plot(ts1,cth125,'-.g','LineWidth',1.25,'color',[0 0.5 0]); hold on
p2 = plot(ts2,cth250,'-b','LineWidth',1.25); hold on
p5 = plot(ts2,cth500,'--r','LineWidth',1.25); hold on
p10 = plot(ts2,cth1000,'--m','LineWidth',1.25); hold on
plot(ts1,(1/1000).*lnb125,'-k','LineWidth',1.25); hold on  % convert from m to km
plot(ts2,(1/1000).*lnb250,'-k','LineWidth',1.25); hold on  % convert from m to km
plot(ts2,(1/1000).*lnb500,'--k','LineWidth',1.25); hold off  % convert from m to km
 
grid on;
box on;
xlabel('Time \big[hr:min\big]','Interpreter','latex');
ylabel('Cloud-top heght \big[km\big]','Interpreter','latex');

% assign legend to specific plots
lgd = legend([p1 p2 p5 p10],{'\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m','\(\Delta=1000\) m'});
set(lgd,'Interpreter','latex','Location','west');

ts_x = 0:5400:length(ts2).*90-1;   % put xticks every 90 mins
xmin = ts_x(1);
xmax = ts_x(end);
xtick = 5400;  % every 90 mins

ymin = 0;
ymax = 12;
ytick = 2;

ax = gca;
ax.XTick = [ts_x(1):xtick:ts_x(end)]; 
ax.XTickLabel = {'00:00','01:30','03:00','04:30','06:00','07:30','09:00','10:30','12:00'};
ax.XLim=[xmin,xmax];
ax.YLim=[ymin,ymax];
ax.YTick = [ymin:ytick:ymax];
% ax.FontSize = 12;
% ax.XLabel.FontSize = 15;
% ax.YLabel.FontSize = 15;
ax.TickLabelInterpreter = 'latex';

pbaspect([2*(xmax-xmin)/3600,ymax-ymin,1]);

% place label at top right corner of the subplot
% text(0.95,1.07,'\big(c\big)','Units','normalized','Interpreter', 'latex','FontWeight','bold'); 

% % place case number on top
% if cn == 1
%     text(0.05,0.9,'CTRL','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end
% if cn == 2
%     text(0.05,0.9,'DRY','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end
% if cn == 3
%     text(0.05,0.9,'MOIST','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end

% END OF CTH PLOT
%-----------------------------------------------------------------------------

% SAVE OUTPUT
fnameg = 'C:\Users\SiuLung\Downloads'; 
filenameg = sprintf('cloudth_1');   % cloud-top height time series

saveas(f,fullfile(fnameg,filenameg),'epsc');
