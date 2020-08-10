close all
clear all

% directory for cloudmf.mat files
% fname = '/aos/home/mtang/Documents/matrices';
% fname = 'C:\Users\SiuLung\Downloads\data_summary';   % for dx = 250
fname = 'C:\Users\SiuLung\Downloads';   % for dx = 250

% enter case number
cn = 1;   % 1 = ctrl, 2 = dry, 3 = moist;

if cn == 1   % ctrl
filename1 = 'ascmf_dx125_may08_wbpqcqi_hydrometeor.mat';
filename2 = 'ascmf_dx250_mar07_wbpqcqi_hydrometeor.mat';
filename5 = 'ascmf_dx500_feb27_wbpqcqi_hydrometeor.mat';
end

if cn == 2   % dry
filename1 = 'ascmf_dx125_dry_jul03_wbpqcqi_hydrometeor.mat';
filename2 = 'ascmf_dx250_dry_jul03_wbpqcqi_hydrometeor.mat';
filename5 = 'ascmf_dx500_dry_jul03_wbpqcqi_hydrometeor.mat';
end

if cn == 3   % moist
filename1 = 'ascmf_dx125_moist_jul03_wbpqcqi_hydrometeor.mat';
filename2 = 'ascmf_dx250_moist_jul03_wbpqcqi_hydrometeor.mat';
filename5 = 'ascmf_dx500_moist_jul03_wbpqcqi_hydrometeor.mat';
end


% directory for output graphs
% fnameg = 'aos/home/mtang/Documents/graphs';
% fnameg = 'C:\Users\SiuLung\Downloads\graphs_summary';
fnameg = 'C:\Users\SiuLung\Downloads';


% % get amf matrcies
amsfx1 = load(fullfile(fname,filename1));  % amsfx array  for dx=125
amsfx2 = load(fullfile(fname,filename2));  % amsfx array  for dx=250
amsfx5 = load(fullfile(fname,filename5));  % amsfx array  for dx=500

fl = ['a' 'b' 'c']; % fig label
P = 1;   % fig position

% control fontsize of texts
lgfs = 11;   % legend font size
axfs = 12.5;   % axis font size
tfs = 13;   % title font size
lbfs = 14;   % x-,y-label font size
cfs = 11;   % case font size (CTRL, DRY, MOIST)
tefs = 15;   % text font size (panel number)

% aspect ratio
pbr = 2.5;

for seq = 6:6

ti = 0+(seq-1)*1; tf = 1+(seq-1)*1;

fieldname = sprintf('f%d%d',ti,tf);

% get z matrix for all plots
z = amsfx2.asmf.(fieldname).z;  % get xh in km

% %----------------------------------------------------------------------------
% % Plot of total ascending mass flux vs z

% get amf from arrays
amf1 = amsfx1.asmf.(fieldname).amf;   % ascending mass flux for dx=125
amf2 = amsfx2.asmf.(fieldname).amf;   % ascending mass flux for dx=250
amf5 = amsfx5.asmf.(fieldname).amf;   % ascending mass flux for dx=500

na1 = amsfx1.asmf.(fieldname).na;   % ascending mass flux for dx=125
na2 = amsfx2.asmf.(fieldname).na;   % ascending mass flux for dx=250
na5 = amsfx5.asmf.(fieldname).na;   % ascending mass flux for dx=500

% % get lcl from arrays
lcl1 = amsfx1.asmf.(fieldname).lclxm;   % lclx for dx=125
lcl2 = amsfx2.asmf.(fieldname).lclxm;   % lclx for dx=250
lcl5 = amsfx5.asmf.(fieldname).lclxm;   % lclx for dx=250

if P == 1
f = figure('units','normalized','outerposition',[0 0 1 1])
end

sp = subplot(1,5,P);

% find mean lclx
lclm = mean(lcl1)/1000;   % in km

% plot(cmf1,lclm*ones(1,length(z)),'--k','LineWidth',1.25); hold on   % mean lcl
plot(amf1,z,'-.g','LineWidth',1.25,'color',[0 0.5 0]); hold on   % dx=125
plot(amf2,z,'-b','LineWidth',1.25); hold on   % dx=250
plot(amf5,z,'--r','LineWidth',1.25); hold on   % dx=500

grid on;
box on;

% xlabel('mass flux \big[kgs\textsuperscript{-1}\big]','Interpreter','latex','FontSize',16);
xlabel('\(M_{c}\) \big[kgs\(^{-1}\)\big]','Interpreter','latex','FontSize',16);

if P == 1
ylabel('Height \big[km\big]','Interpreter','latex','FontSize',16);
end

xmin = 0;
xmax = 3e8;
xtick = 1e8;

plot([xmin xmax],[lclm lclm],'--k','LineWidth',1.25); hold off

% lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m','LCL');
% set(lgd,'Interpreter','latex','Location','northeast');
% lgd.FontSize = lgfs;   % 14 originally

ymin = 0;
ymax = 10;  
ytick = 2;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];
ax.XAxis.Exponent = 9;

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];
% ax.YAxis.Exponent = 4;
ax.TickLabelInterpreter = 'latex';

ax.FontSize = axfs;  % 16 originally
ax.XLabel.FontSize = lbfs;
ax.YLabel.FontSize = lbfs;

% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),pbr*(xmax-xmin),1]); % multiple y-axis by the factor

tih = ti; timin = 0;
tfh = tf; tfmin = 0;

timehi = sprintf('%02d',tih);
timemini = sprintf('%02d',timin);
timehf = sprintf('%02d',tfh);
timeminf = sprintf('%02d',tfmin);

% title([num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf) ':' ...
%     num2str(timeminf)],'Interpreter','latex','FontSize',tfs);


% place labels on plot
text(0.88,1.03,'(a)','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

if cn == 1
    text(0.1,0.95,'CTRL','Units','normalized','Interpreter', 'latex','FontSize',cfs,'FontWeight','bold');
end
if cn == 2
    text(0.1,0.95,'DRY','Units','normalized','Interpreter', 'latex','FontSize',cfs,'FontWeight','bold');
end
if cn == 3
    text(0.1,0.95,'MOIST','Units','normalized','Interpreter', 'latex','FontSize',cfs,'FontWeight','bold');
end

% save axes info into varaible for editing space between subplots later
ax1 = ax;

% END OF MASS FLUX OF ASCENDING POINTS PLOT
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

% PLOT of profile of total area of ascending points

% get w_am from arrays
aat1 = amsfx1.asmf.(fieldname).aat;   % total area in m^2 of asc. points for dx=125
aat2 = amsfx2.asmf.(fieldname).aat;   % total area in m^2 of asc. points for dx=250
aat5 = amsfx5.asmf.(fieldname).aat;   % total area in m^2 of asc. points for dx=500

% convert aat to km^2
aat1k = aat1./(1000^2);  % in km^2
aat2k = aat2./(1000^2);  % in km^2
aat5k = aat5./(1000^2);  % in km^2

sp = subplot(1,5,2);

plot(aat1k,z,'-.g','LineWidth',1.25,'color',[0 0.5 0]); hold on   % dx=125
plot(aat2k,z,'-b','LineWidth',1.25); hold on   % dx=250
plot(aat5k,z,'--r','LineWidth',1.25); hold on   % dx=500

% plot([xmin xmax],[lclm lclm],'--k','LineWidth',1.25); hold off

grid on;
box on;

xlabel('\(A\) \big[km\(^{2}\)]','Interpreter','latex','FontSize',16);
% xlabel('total area \big[km\textsuperscript{2}]','Interpreter','latex','FontSize',16);

xmin = 0;
xmax = 150;
xtick = 50;

plot([xmin xmax],[lclm lclm],'--k','LineWidth',1.25); hold off

ymin = 0;
ymax = 10;  
ytick = 2;

% lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m','LCL');
% set(lgd,'Interpreter','latex','Location','northeast');
% lgd.FontSize = lgfs;   % 14 originally

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];
% ax.XAxis.Exponent = 3;

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];
% ax.YAxis.Exponent = 4;
ax.TickLabelInterpreter = 'latex';

ax.FontSize = axfs;  % 16 originally
ax.XLabel.FontSize = lbfs;

% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),pbr*(xmax-xmin),1]); % multiple y-axis by the factor

tih = ti; timin = 0;
tfh = tf; tfmin = 0;

timehi = sprintf('%02d',tih);
timemini = sprintf('%02d',timin);
timehf = sprintf('%02d',tfh);
timeminf = sprintf('%02d',tfmin);

% title([num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf) ':' ...
%     num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% place labels on plot
text(0.88,1.03,'(b)','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% if cn == 1
%     text(0.7,0.60,'CTRL','Units','normalized','Interpreter', 'latex','FontSize',cfs,'FontWeight','bold');
% end
% if cn == 2
%     text(0.1,0.95,'DRY','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end
% if cn == 3
%     text(0.1,0.95,'MOIST','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end

% text(0.70,0.55,'\(w>0\)','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(0.70,0.50,'\(b>0\)','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(0.38,0.45,'\(q_c+q_i>0.1 \) g/kg','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% 
% text(0.35,0.40,'\(-60\leq x\leq 60\) km','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% save axes info into varaible for editing space between subplots later
ax2 = ax;

% END OF TOTAL AREA OF ASCENDING POINTS PLOT
%----------------------------------------------------------------------------

% PLOT of avg. vertical velocity vs z

% get w_am from arrays
wam1 = amsfx1.asmf.(fieldname).wam;   % avg. vertical velocity for dx=125
wam2 = amsfx2.asmf.(fieldname).wam;   % avg. vertical velocity for dx=250
wam5 = amsfx5.asmf.(fieldname).wam;   % avg. vertical velocity for dx=500

sp = subplot(1,5,3);

plot(wam1,z,'-.g','LineWidth',1.25,'color',[0 0.5 0]); hold on   % dx=125
plot(wam2,z,'-b','LineWidth',1.25); hold on   % dx=250
plot(wam5,z,'--r','LineWidth',1.25); hold on   % dx=500

plot([xmin xmax],[lclm lclm],'--k','LineWidth',1.25); hold off

grid on;
box on;

% lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m','LCL');
% set(lgd,'Interpreter','latex','Location','northeast');
% lgd.FontSize = 12;   % 14 originally

xlabel('\(w_c\) \big[ms\(^{-1}\)]','Interpreter','latex','FontSize',16);
% xlabel('\(\overline{w}\) \big[ms\textsuperscript{-1}]','Interpreter','latex','FontSize',16);

xmin = 0;
xmax = 10;
xtick = 2;

ymin = 0;
ymax = 10;  
ytick = 2;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];
% ax.YAxis.Exponent = 4;
ax.TickLabelInterpreter = 'latex';

ax.FontSize = axfs;  % 16 originally
ax.XLabel.FontSize = lbfs;

% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),pbr*(xmax-xmin),1]); % multiple y-axis by the factor

tih = ti; timin = 0;
tfh = tf; tfmin = 0;

timehi = sprintf('%02d',tih);
timemini = sprintf('%02d',timin);
timehf = sprintf('%02d',tfh);
timeminf = sprintf('%02d',tfmin);

% title([num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf) ':' ...
%     num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% place labels on plot
text(0.88,1.03,'(c)','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% if cn == 1
%     text(0.1,0.95,'CTRL','Units','normalized','Interpreter', 'latex','FontSize',cfs,'FontWeight','bold');
% end
% if cn == 2
%     text(0.1,0.95,'DRY','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end
% if cn == 3
%     text(0.1,0.95,'MOIST','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end

% save axes info into varaible for editing space between subplots later
ax3 = ax;

% END OF AVG. VERTICAL VELOCITY OF ASCENDING POINTS PLOT

%----------------------------------------------------------------------------

% PLOT of profile of avg. buoyancy of ascending points

% get w_am from arrays
bam1 = amsfx1.asmf.(fieldname).bam;   % avg. buoyancy of asc. points for dx=125
bam2 = amsfx2.asmf.(fieldname).bam;   % avg. buoyancy of asc. points for dx=250
bam5 = amsfx5.asmf.(fieldname).bam;   % avg. buoyancy of asc. points for dx=500

sp = subplot(1,5,4);

plot(bam1,z,'-.g','LineWidth',1.25,'color',[0 0.5 0]); hold on   % dx=125
plot(bam2,z,'-b','LineWidth',1.25); hold on   % dx=250
plot(bam5,z,'--r','LineWidth',1.25); hold on   % dx=500

% plot([xmin xmax],[lclm lclm],'--k','LineWidth',1.25); hold off

grid on;
box on;

xlabel('\(b_c\) \big[ms\(^{-2}\)]','Interpreter','latex','FontSize',16);
% xlabel('\(\overline{b}\) \big[ms\textsuperscript{-2}]','Interpreter','latex','FontSize',16);

xmin = 0;
xmax = 2.5e-2;
xtick = 0.5e-2; 

plot([xmin xmax],[lclm lclm],'--k','LineWidth',1.25); hold off

ymin = 0;
ymax = 10;  
ytick = 2;

% lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m','LCL');
% set(lgd,'Interpreter','latex','Location','northeast');
% lgd.FontSize = 12;   % 14 originally

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];
ax.XAxis.Exponent = -2;

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];
% ax.YAxis.Exponent = 4;
ax.TickLabelInterpreter = 'latex';

ax.FontSize = axfs;  % 16 originally
ax.XLabel.FontSize = lbfs;

% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),pbr*(xmax-xmin),1]); % multiple y-axis by the factor

tih = ti; timin = 0;
tfh = tf; tfmin = 0;

timehi = sprintf('%02d',tih);
timemini = sprintf('%02d',timin);
timehf = sprintf('%02d',tfh);
timeminf = sprintf('%02d',tfmin);

% title([num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf) ':' ...
%     num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% place labels on plot
text(0.88,1.03,'(d)','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% if cn == 1
%     text(0.1,0.95,'CTRL','Units','normalized','Interpreter', 'latex','FontSize',cfs,'FontWeight','bold');
% end
% if cn == 2
%     text(0.1,0.95,'DRY','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end
% if cn == 3
%     text(0.1,0.95,'MOIST','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end
% 
% text(0.05,0.45,'\(w>0\)','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(0.05,0.40,'\(b>0\)','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(0.05,0.35,'\(q_c+q_i>0.1 \) g/kg','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% save axes info into varaible for editing space between subplots later
ax4 = ax;

% END OF AVG. BUOYANCY OF ASCENDING POINTS PLOT
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

% PLOT of profile of avg. qh of ascending points

% get qt_am from arrays
qham1 = amsfx1.asmf.(fieldname).qham;   % in kg/kg, avg. qh of asc. points for dx=125
qham2 = amsfx2.asmf.(fieldname).qham;   % in kg/kg, avg. qh of asc. points for dx=250
qham5 = amsfx5.asmf.(fieldname).qham;   % in kg/kg, avg. qh of asc. points for dx=500

sp = subplot(1,5,5);

% convert qh from kg/kg to g/kg
qham1g = 1000.*qham1;
qham2g = 1000.*qham2;
qham5g = 1000.*qham5;

plot(qham1g,z,'-.g','LineWidth',1.25,'color',[0 0.5 0]); hold on   % dx=125
plot(qham2g,z,'-b','LineWidth',1.25); hold on   % dx=250
plot(qham5g,z,'--r','LineWidth',1.25); hold on   % dx=500

% plot([xmin xmax],[lclm lclm],'--k','LineWidth',1.25); hold off

grid on;
box on;

% xlabel('\(q_{hc}\) \big[gkg\textsuperscript{-1}]','Interpreter','latex','FontSize',16);
xlabel('\(q_{hc}\) \big[gkg\(^{-1}\)]','Interpreter','latex','FontSize',16);
% xlabel('\(\overline{q_h}\) \big[gkg\textsuperscript{-1}]','Interpreter','latex','FontSize',16);

xmin = 0;
xmax = 3;
xtick = 1; 

plot([xmin xmax],[lclm lclm],'--k','LineWidth',1.25); hold off

lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m','LCL');
set(lgd,'Interpreter','latex','Location','northeast');
lgd.FontSize = lgfs;   % 14 originally

ymin = 0;
ymax = 10;  
ytick = 2;

% lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m','LCL');
% set(lgd,'Interpreter','latex','Location','northeast');
% lgd.FontSize = 14;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];
% ax.XAxis.Exponent = -2;

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];
% ax.YAxis.Exponent = 4;
ax.TickLabelInterpreter = 'latex';

ax.FontSize = axfs;   % 26 originally
ax.XLabel.FontSize = lbfs;

% % Make y-axis:x-axis = 2
pbaspect([(xmax-xmin),pbr*(xmax-xmin),1]); % multiple y-axis by the factor

tih = ti; timin = 0;
tfh = tf; tfmin = 0;

timehi = sprintf('%02d',tih);
timemini = sprintf('%02d',timin);
timehf = sprintf('%02d',tfh);
timeminf = sprintf('%02d',tfmin);

% title([num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf) ':' ...
%     num2str(timeminf)],'Interpreter','latex','FontSize',tfs);

% place labels on plot
text(0.88,1.03,'(e)','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% if cn == 1
%     text(0.1,0.95,'CTRL','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end
% if cn == 2
%     text(0.1,0.95,'DRY','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end
% if cn == 3
%     text(0.1,0.95,'MOIST','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end

% if cn == 1
%     text(0.7,0.60,'CTRL','Units','normalized','Interpreter', 'latex','FontSize',cfs,'FontWeight','bold');
% end
% if cn == 2
%     text(0.75,0.50,'DRY','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end
% if cn == 3
%     text(0.75,0.50,'MOIST','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end

% text(0.55,0.45,'\(w>0\)','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');
% text(0.35,0.40,'\(q_c+q_i>0.1 \) g/kg','Units','normalized','Interpreter', 'latex','FontSize',tefs,'FontWeight','bold');

% save axes info into varaible for editing space between subplots later
ax5 = ax;

% END OF AVG. QT PLOT
%----------------------------------------------------------------------------

% PLOT of total number of ascending points

% % get total number of ascending points from arrays
% na1 = amsfx1.asmf.(fieldname).na;  % dx=125
% na2 = amsfx2.asmf.(fieldname).na;  % dx=250
% na5 = amsfx5.asmf.(fieldname).na;  % dx=500
% 
% figure
% plot(na1,z,'-.g','LineWidth',1.25,'color',[0 0.5 0]); hold on   % dx=125
% plot(na2,z,'-b','LineWidth',1.25); hold on   % dx=250
% plot(na5,z,'--r','LineWidth',1.25); hold off   % dx=500
% % 



% % END OF cloud-bse mass flux vs x plot
% %----------------------------------------------------------------------------
% %----------------------------------------------------------------------------

% Reduce space between subplots

ax2.Position(3) = ax1.Position(3);
ax3.Position(3) = ax1.Position(3);
ax4.Position(3) = ax1.Position(3);
ax5.Position(3) = ax1.Position(3);

ax2.Position(4) = ax1.Position(4);
ax3.Position(4) = ax1.Position(4);
ax4.Position(4) = ax1.Position(4);
ax5.Position(4) = ax1.Position(4);

ax2.Position(1) = ax1.Position(1) + ax1.Position(3) + 0.025;
ax3.Position(1) = ax2.Position(1) + ax2.Position(3) + 0.025;
ax4.Position(1) = ax3.Position(1) + ax3.Position(3) + 0.025;
ax5.Position(1) = ax4.Position(1) + ax4.Position(3) + 0.025;

ax1.Position(3) = ax2.Position(3);



% SAVE FIGURE

if cn == 1
filenameg = sprintf('ascmf_wbpqcqi_%d%d_ctrl_f',ti,tf);
end

if cn == 2
filenameg = sprintf('ascmf_wbpqcqi_%d%d_dry',ti,tf);
end

if cn == 3
filenameg = sprintf('ascmf_wbpqcqi_%d%d_moist',ti,tf);
end

% saveas(sp,fullfile(fnameg,filenameg),'epsc');


end   % end seq = 6:6
