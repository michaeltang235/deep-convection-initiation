close all
clear all

tic 

% directory for input files
% fname = '/aos/home/mtang/Documents/matrices';
fname = 'C:\Users\SiuLung\Downloads\paper_data';

% enter case number
cn = 3; % 1 = ctrl, 2 = dry, 3 = moist

% get the corresponding filenames
if cn == 1
filename1 = 'convm_dx125_may08_2.txt';
filename2 = 'convm_dx250_mar07_2.txt';
filename5 = 'convm_dx500_feb27_2.txt';
end

if cn == 2
filename1 = 'convm_dx125_dry_jul03_2.txt';
filename2 = 'convm_dx250_dry_jul03_2.txt';
filename5 = 'convm_dx500_dry_jul03_2.txt';
end

if cn == 3
filename1 = 'convm_dx125_moist_jul03_2.txt';
filename2 = 'convm_dx250_moist_jul03_2.txt';
filename5 = 'convm_dx500_moist_jul03_2.txt';
end


convm125 = load(fullfile(fname,filename1));  % dx = 125
convm250 = load(fullfile(fname,filename2));  % dx = 250
convm500 = load(fullfile(fname,filename5));  % dx = 500

% convm250 = dlmread('/aos/home/mtang/Documents/matrices/convm_dx250_mar07_2.txt');  % dx = 250
% convm500 = dlmread('/aos/home/mtang/Documents/matrices/convm_dx500_feb27_2.txt');  % dx = 500

ts1 = ((1:1:121)-1).*360;  % make time series for 121 input files
ts2 = ((1:1:481)-1).*90;  % make time series for 481 input files

%--------------------------------------------------------
% CONVM PLOT
% f1 = figure
% open figure window at screen size and at normalized units
f = figure('units','normalized','outerposition',[0 0 1 1])
% set(gca,'FontSize',9,'LineWidth',1.5);
sp = subplot(2,2,1)
p1 = plot(ts1,convm125,'-.g','LineWidth',1.25,'color',[0 0.5 0]); hold on
p2 = plot(ts2,convm250,'-b','LineWidth',1.25); hold on
p5 = plot(ts2,convm500,'--r','LineWidth',1.25); hold off

grid on;
box on;
xlabel('Time \big[hr:min\big]','Interpreter','latex','FontSize',16);
ylabel('Convergence \big[\(s^{-1}\)\big]','Interpreter','latex','FontSize',16);

lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m');
set(lgd,'Interpreter','latex','Location','southeast');


ts_x = 0:5400:length(ts2).*90-1;   % put xticks every 90 mins
xmin = ts_x(1);
xmax = ts_x(end);
xlim([xmin xmax]);
% xtick = 3600;
xtick = 5400;  % every 90 mins

ymin = -2/1000;
ymax = 10/1000;  
ytick = 2/1000;

ax = gca;
ax.XTick = [ts_x(1):xtick:ts_x(end)]; 
ax.XTickLabel = {'00:00','01:30','03:00','04:30','06:00','07:30','09:00','10:30','12:00'};
ax.YLim=[ymin,ymax];
ax.YTick = [ymin:ytick:ymax];
ax.YAxis.Exponent = -3;
ax.FontSize = 12;
ax.XLabel.FontSize = 15;
ax.YLabel.FontSize = 15;
ax.TickLabelInterpreter = 'latex';
pbaspect([2*(xmax-xmin)/3600,1000*(ymax-ymin),1]);

% place label at top right corner of the subplot
text(0.95,1.07,'\big(a\big)','Units','normalized','Interpreter', 'latex','FontSize',16,'FontWeight','bold'); 

if cn == 1
    text(0.05,0.9,'CTRL','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
end
if cn == 2
    text(0.05,0.9,'DRY','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
end
if cn == 3
    text(0.05,0.9,'MOIST','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
end

% assign current axes object ax1 (for adjusting spacing btw subplots later)
ax1 = ax;

% END OF CONVM PLOT
%--------------------------------------------------------
%
%--------------------------------------------------------
% CAPE_CIN_PLOT

if cn == 1
filenamecp1 = 'capem_dx125_may08.txt';
filenamecp2 = 'capem_dx250_mar07.txt';
filenamecp5 = 'capem_dx500_feb27.txt';

filenameci1 = 'cinm_dx125_may08.txt';
filenameci2 = 'cinm_dx250_mar07.txt';
filenameci5 = 'cinm_dx500_feb27.txt';
end

if cn == 2
filenamecp1 = 'capem_dx125_dry_jul03.txt';
filenamecp2 = 'capem_dx250_dry_jul03.txt';
filenamecp5 = 'capem_dx500_dry_jul03.txt';

filenameci1 = 'cinm_dx125_dry_jul03.txt';
filenameci2 = 'cinm_dx250_dry_jul03.txt';
filenameci5 = 'cinm_dx500_dry_jul03.txt';
end

if cn == 3
filenamecp1 = 'capem_dx125_moist_jul03.txt';
filenamecp2 = 'capem_dx250_moist_jul03.txt';
filenamecp5 = 'capem_dx500_moist_jul03.txt';

filenameci1 = 'cinm_dx125_moist_jul03.txt';
filenameci2 = 'cinm_dx250_moist_jul03.txt';
filenameci5 = 'cinm_dx500_moist_jul03.txt';
end


capem125 = load(fullfile(fname,filenamecp1));  % dx = 125
cinm125 = load(fullfile(fname,filenameci1));  % dx = 125

capem250 = load(fullfile(fname,filenamecp2));  % dx = 250
cinm250 = load(fullfile(fname,filenameci2));  % dx = 250

capem500 = load(fullfile(fname,filenamecp5));  % dx = 500
cinm500 = load(fullfile(fname,filenameci5));  % dx = 500

% smoothen curve by removing anormalies 
% taking avg. of adjacnet entries
tol = 10;
cinm125 = smoothing(cinm125,tol);
cinm250 = smoothing(cinm250,tol);
cinm500 = smoothing(cinm500,tol);

ts2 = ((1:1:481)-1).*90;  % make time series for 481 input files

sp = subplot(2,2,2);
[hAx3, p5, p6]=plotyy(ts1,capem125,ts1,cinm125); hold on
[hAx1, p1, p2]=plotyy(ts2,capem250,ts2,cinm250); hold on
[hAx2, p3, p4]=plotyy(ts2,capem500,ts2,cinm500); hold off

grid on;
box on;

% change line style for each plot using handles 
p1.LineStyle = '-';   % dx=250
p2.LineStyle = '-';
p3.LineStyle = '--';   % dx=500
p4.LineStyle = '--';
p5.LineStyle = '-.';   % dx=125
p6.LineStyle = '-.';

% p1.Marker = '+';

% change line width for each plot using handles 
p1.LineWidth = 1.25;
p2.LineWidth = 1.25;
p3.LineWidth = 1.25;
p4.LineWidth = 1.25;
p5.LineWidth = 1.25;
p6.LineWidth = 1.25;

% change line color for each plot using handles 
p1.Color = 'b';
p2.Color = 'b';
p3.Color = 'r';
p4.Color = 'r';
p5.Color = [0 0.5 0];
p6.Color = [0 0.5 0];

% assign legend to specific plots
% legend([p1 p3],{'dx=250','dx=500'});
lgd = legend([p5 p1 p3],{'\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m'});
set(lgd,'Interpreter','latex','Location','northeast');

% format x-axis
ts_x = 0:5400:length(ts2).*90-1;   % put xticks every 90 mins
xmin = ts_x(1);
xmax = ts_x(end);
xlim([xmin xmax]);
xtick = 5400;  % every 90 mins

% set x-axis limit
hAx1(1).XLim = [xmin,xmax];
hAx1(2).XLim = [xmin,xmax];
hAx2(1).XLim = [xmin,xmax];
hAx2(2).XLim = [xmin,xmax];
hAx3(1).XLim = [xmin,xmax];
hAx3(2).XLim = [xmin,xmax];

hAx1(1).XTick = [ts_x(1):xtick:ts_x(end)]; 
hAx1(1).XTickLabel = {'00:00','01:30','03:00','04:30','06:00','07:30','09:00','10:30','12:00'};

hAx1(2).XTick = [ts_x(1):xtick:ts_x(end)]; 
hAx1(2).XTickLabel = {'00:00','01:30','03:00','04:30','06:00','07:30','09:00','10:30','12:00'};

hAx2(1).XTick = [ts_x(1):xtick:ts_x(end)]; 
hAx2(1).XTickLabel = {'00:00','01:30','03:00','04:30','06:00','07:30','09:00','10:30','12:00'};

hAx2(2).XTick = [ts_x(1):xtick:ts_x(end)]; 
hAx2(2).XTickLabel = {'00:00','01:30','03:00','04:30','06:00','07:30','09:00','10:30','12:00'};

hAx3(1).XTick = [ts_x(1):xtick:ts_x(end)]; 
hAx3(1).XTickLabel = {'00:00','01:30','03:00','04:30','06:00','07:30','09:00','10:30','12:00'};

% format both y-axes
ymin = 900;  % range in left y-axis
ymax = 1700;  
ytick = 200;

ymin2 = 0;  % range in right y-axis
ymax2 = 60;  
ytick2 = 20;

hAx1(1).YLim=[ymin,ymax];   % set limit of left y-axis
hAx1(1).YTick = [ymin:ytick:ymax];

hAx2(1).YLim=[ymin,ymax];   % set limit of left y-axis
hAx2(1).YTick = [ymin:ytick:ymax];

hAx3(1).YLim=[ymin,ymax];   % set limit of left y-axis
hAx3(1).YTick = [ymin:ytick:ymax];


hAx1(2).YLim=[ymin2,ymax2];   % set limit of right y-axis
hAx1(2).YTick = [ymin2:ytick2:ymax2];

hAx2(2).YLim=[ymin2,ymax2];   % set limit of right y-axis
hAx2(2).YTick = [ymin2:ytick2:ymax2];

hAx3(2).YLim=[ymin2,ymax2];   % set limit of right y-axis
hAx3(2).YTick = [ymin2:ytick2:ymax2];

hAx1(1).FontSize = 12;
hAx1(2).FontSize = 12;
hAx2(1).FontSize = 12;
hAx2(2).FontSize = 12;
hAx3(1).FontSize = 12;
hAx3(2).FontSize = 12;


% set both y-axes to black color
hAx1(1).YColor = 'k';
hAx1(2).YColor = 'k';
hAx2(1).YColor = 'k';
hAx2(2).YColor = 'k';
hAx3(1).YColor = 'k';
hAx3(2).YColor = 'k';

% set axes interpreter as latex
hAx1(1).TickLabelInterpreter = 'latex';


% set ratio of x-axis:y-axis = 2
pbaspect(hAx1(1),[2*(ymax-ymin),(ymax-ymin),1]);
pbaspect(hAx1(2),[2*(ymax2-ymin2),(ymax2-ymin2),1]);
pbaspect(hAx2(1),[2*(ymax-ymin),(ymax-ymin),1]);
pbaspect(hAx2(2),[2*(ymax2-ymin2),(ymax2-ymin2),1]);
pbaspect(hAx3(1),[2*(ymax-ymin),(ymax-ymin),1]);
pbaspect(hAx3(2),[2*(ymax2-ymin2),(ymax2-ymin2),1]);

% format x- and y- labels
xlabel('Time \big[hr:min\big]','Interpreter','latex','FontSize',15);
ylabel(hAx1(1),'CAPE \big[\(Jkg^{-1}\)\big]','Interpreter','latex','FontSize',15); % left y-axis 
ylabel(hAx1(2),'CIN \big[\(Jkg^{-1}\)\big]','Interpreter','latex','FontSize',15); % right y-axis

% place label at top right corner of the subplot
text(0.95,1.07,'\big(b\big)','Units','normalized','Interpreter', 'latex','FontSize',16,'FontWeight','bold'); 

% place case number on top
if cn == 1
    text(0.05,0.9,'CTRL','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
end
if cn == 2
    text(0.05,0.9,'DRY','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
end
if cn == 3
    text(0.05,0.9,'MOIST','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
end

% END OF CAPE_CIN_PLOT
%----------------------------------------------------------------------------

%----------------------------------------------------------------------------
% PLOT CLOUD TOP HEIGHT 

if cn == 1
filenamect1 = 'cth_dx125_may08_1.txt';
filenamel1 = 'lnb_dx125_may08.txt';

filenamect2 = 'cth_dx250_mar07_1.txt';
filenamel2 = 'lnb_dx250_mar07.txt';

filenamect5 = 'cth_dx500_feb27_1.txt';
filenamel5 = 'lnb_dx500_feb27.txt';
end

if cn == 2
filenamect1 = 'cth_dx125_dry_jul03_1.txt';
filenamel1 = 'lnb_dx125_dry_jul03.txt';

filenamect2 = 'cth_dx250_dry_jul03_1.txt';
filenamel2 = 'lnb_dx250_dry_jul03.txt';

filenamect5 = 'cth_dx500_dry_jul03_1.txt';
filenamel5 = 'lnb_dx500_dry_jul03.txt';
end

if cn == 3
filenamect1 = 'cth_dx125_moist_jul03_1.txt';
filenamel1 = 'lnb_dx125_moist_jul03.txt';

filenamect2 = 'cth_dx250_moist_jul03_1.txt';
filenamel2 = 'lnb_dx250_moist_jul03.txt';

filenamect5 = 'cth_dx500_moist_jul03_1.txt';
filenamel5 = 'lnb_dx500_moist_jul03.txt';
end

cth125 = load(fullfile(fname,filenamect1));
lnb125 = load(fullfile(fname,filenamel1));

cth250 = load(fullfile(fname,filenamect2));
lnb250 = load(fullfile(fname,filenamel2));

cth500 = load(fullfile(fname,filenamect5));
lnb500 = load(fullfile(fname,filenamel5));

% % get avg. lnb
% lnb125 = mean(mean(lnb125));
% lnb250 = mean(mean(lnb250));
% lnb500 = mean(mean(lnb500));

% cth500 = smoothing(cth500,3);

sp = subplot(2,2,3);

p1 = plot(ts1,cth125,'-.g','LineWidth',1.25,'color',[0 0.5 0]); hold on
p2 = plot(ts2,cth250,'-b','LineWidth',1.25); hold on
p5 = plot(ts2,cth500,'--r','LineWidth',1.25); hold on
plot(ts1,(1/1000).*lnb125,'-k','LineWidth',1.25); hold on  % convert from m to km
plot(ts2,(1/1000).*lnb250,'-k','LineWidth',1.25); hold on  % convert from m to km
plot(ts2,(1/1000).*lnb500,'--k','LineWidth',1.25); hold off  % convert from m to km
 
grid on;
box on;
xlabel('Time \big[hr:min\big]','Interpreter','latex','FontSize',16);
ylabel('Cloud-top heght \big[km\big]','Interpreter','latex','FontSize',16);

% assign legend to specific plots
lgd = legend([p1 p2 p5],{'\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m'});
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
ax.FontSize = 12;
ax.XLabel.FontSize = 15;
ax.YLabel.FontSize = 15;
ax.TickLabelInterpreter = 'latex';

pbaspect([2*(xmax-xmin)/3600,ymax-ymin,1]);

% place label at top right corner of the subplot
text(0.95,1.07,'\big(c\big)','Units','normalized','Interpreter', 'latex','FontSize',16,'FontWeight','bold'); 

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

% assign current axes object ax3 (for adjusting spacing btw subplots later)
ax3 = ax;

% END OF CTH PLOT
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% PLOT RELATIVE HUMIDITY

if cn == 1
filenamer1 = 'rhm_dx125_may08_3.txt';
filenamer2 = 'rhm_dx250_mar07_3.txt';
filenamer5 = 'rhm_dx500_feb27_3.txt';

% filenamer1 = 'rhm_dx125_may08_2.txt';
% filenamer2 = 'rhm_dx250_mar07_2.txt';
% filenamer5 = 'rhm_dx500_feb27_2.txt';
end

if cn == 2
filenamer1 = 'rhm_dx125_dry_jul03_3.txt';
filenamer2 = 'rhm_dx250_dry_jul03_3.txt';
filenamer5 = 'rhm_dx500_dry_jul03_3.txt';

% filenamer1 = 'rhm_dx125_dry_jul03_2.txt';
% filenamer2 = 'rhm_dx250_dry_jul03_2.txt';
% filenamer5 = 'rhm_dx500_dry_jul03_2.txt';
end

if cn == 3
filenamer1 = 'rhm_dx125_moist_jul03_3.txt';
filenamer2 = 'rhm_dx250_moist_jul03_3.txt';
filenamer5 = 'rhm_dx500_moist_jul03_3.txt';
% filenamer1 = 'rhm_dx125_moist_jul03_2.txt';
% filenamer2 = 'rhm_dx250_moist_jul03_2.txt';
% filenamer5 = 'rhm_dx500_moist_jul03_2.txt';
end

rhm125 = load(fullfile(fname,filenamer1));  % dx = 125
rhm250 = load(fullfile(fname,filenamer2));  % dx = 250
rhm500 = load(fullfile(fname,filenamer5));  % dx = 250

sp = subplot(2,2,4);

plot(ts1,rhm125,'-.g','LineWidth',1.25,'color',[0 0.5 0]); hold on
plot(ts2,rhm250,'-b','LineWidth',1.25); hold on
plot(ts2,rhm500,'--r','LineWidth',1.25); hold off

grid on;
box on;
xlabel('Time \big[hr:min\big]','Interpreter','latex','FontSize',16);
ylabel('RH mid-troposphere','Interpreter','latex','FontSize',16);

lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m');
if cn == 1 || cn == 3
set(lgd,'Interpreter','latex','Location','southeast');
else 
    set(lgd,'Interpreter','latex','Location','northeast');
end
% lgd.FontSize = 14;

% ts_x = 0:3600:length(ts2).*90-1;  % put xticks every 60 mins
ts_x = 0:5400:length(ts2).*90-1;   % put xticks every 90 mins
xmin = ts_x(1);
xmax = ts_x(end);
xlim([xmin xmax]);
xtick = 5400;  % every 90 mins

% ymin = 0.67; % originally is 0.5
% ymax = 0.71; % originally is 0.9 
% ytick = 0.01;  % originally is 0.1

ymin = 0.65; % originally is 0.5
ymax = 0.75; % originally is 0.9 
ytick = 0.025;  % originally is 0.1

if cn == 2
    ymin = 0.55; % originally is 0.5
    ymax = 0.65; % originally is 0.9 
    ytick = 0.025;  % originally is 0.1
end

if cn == 3
    ymin = 0.75; % originally is 0.5
    ymax = 0.85; % originally is 0.9 
    ytick = 0.025;  % originally is 0.1
end

% if cn == 2 || cn == 3
%     ymin = 0.55; % originally is 0.5
%     ymax = 0.85; % originally is 0.9 
%     ytick = 0.1;  % originally is 0.1
% end

ax = gca;
ax.XTick = [ts_x(1):xtick:ts_x(end)]; 
% ax.XTickLabel = {'00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00','08:00','09:00','10:00','11:00','12:00'};
ax.XTickLabel = {'00:00','01:30','03:00','04:30','06:00','07:30','09:00','10:30','12:00'};
ax.YLim=[ymin,ymax];
ax.YTick = [ymin:ytick:ymax];
ax.FontSize = 12;
ax.XLabel.FontSize = 15;
ax.YLabel.FontSize = 15;
ax.TickLabelInterpreter = 'latex';

% make x-axis:y-axis = 2:1
pbaspect([2*(ymax-ymin),ymax-ymin,1]);

% place label at top right corner of the subplot
text(0.95,1.07,'\big(d\big)','Units','normalized','Interpreter', 'latex','FontSize',16,'FontWeight','bold'); 

% place case number on top
if cn == 1
    text(0.05,0.9,'CTRL','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
end
if cn == 2
    text(0.05,0.9,'DRY','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
end
if cn == 3
    text(0.05,0.9,'MOIST','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
end

% assign current axes object ax3 (for adjusting spacing btw subplots later)
ax4 = ax;

% END OF RH PLOT
%----------------------------------------------------------------------------

% reduce spacing between subplots
ax1.Position(1) = hAx1(1).Position(1) - ax1.Position(3)-0.065;
ax3.Position(1) = ax1.Position(1);

ax3.Position(2) = ax1.Position(2) - ax1.Position(4) - 0.095;
ax4.Position(2) = ax3.Position(2);

%----------------------------------------------------------------------------

% SAVE FIGURE
% fname = '/aos/home/mtang/Documents/graphs';
fnameg = 'C:\Users\SiuLung\Downloads';

if cn == 1
filenameg = 'subplot_4fig_1_ctrl';
end

if cn == 2
filenameg = 'subplot_4fig_1_dry';
end

if cn == 3
filenameg = 'subplot_4fig_1_moist';
end

% saveas(f,fullfile(fnameg,filenameg),'epsc');

% saveas(f,'subplot_4fig','epsc')
% print(gcf,'subplot_4fig','-depsc','-r0')

toc
