close all
clear all

%----------------------------------------------------------------------------
% This code generates plot of frac. ent. and dent. rates from data
% calculated using ice-iquid water potential temp.
%----------------------------------------------------------------------------

% directory for msfx,fmsfx,entm, and detm .mat files
% fname = 'aos/home/mtang/Documents/matrices';
fname = 'C:\Users\SiuLung\Downloads';

% enter case number
cn = 1;   % 1 = ctrl, 2 = dry, 3 = moist;

% get corresponding filenames

if cn == 1
filename01 = 'msfx_dx125_may08_4.mat';   % for dx=125, ctrl
filename21 = 'entm_thil_dx125_may08.mat';
filename31 = 'detm_thil_dx125_may08.mat';

filename0 = 'msfx_dx250_mar07_4.mat';   % for dx=250, ctrl
filename2 = 'entm_thil_dx250_mar07.mat';
filename3 = 'detm_thil_dx250_mar07.mat';

filename05 = 'msfx_dx500_feb27_4.mat';   % for dx = 500, ctrl
filename25 = 'entm_thil_dx500_feb27.mat';
filename35 = 'detm_thil_dx500_feb27.mat';
end

if cn == 2
filename01 = 'msfx_dx125_dry_jul03_4.mat';   % for dx=125, dry
filename21 = 'entm_thil_dx125_dry_jul03.mat';
filename31 = 'detm_thil_dx125_dry_jul03.mat';

filename0 = 'msfx_dx250_dry_jul03_4.mat';   % for dx=250, dry
filename2 = 'entm_thil_dx250_dry_jul03.mat';
filename3 = 'detm_thil_dx250_dry_jul03.mat';

filename05 = 'msfx_dx500_dry_jul03_4.mat';   % for dx = 500, dry
filename25 = 'entm_thil_dx500_dry_jul03.mat';
filename35 = 'detm_thil_dx500_dry_jul03.mat';
end

if cn == 3
filename01 = 'msfx_dx125_moist_jul03_4.mat';   % for dx=125, moist
filename21 = 'entm_thil_dx125_moist_jul03.mat';
filename31 = 'detm_thil_dx125_moist_jul03.mat';

filename0 = 'msfx_dx250_moist_jul03_4.mat';   % for dx=250, moist
filename2 = 'entm_thil_dx250_moist_jul03.mat';
filename3 = 'detm_thil_dx250_moist_jul03.mat';

filename05 = 'msfx_dx500_moist_jul03_4.mat';   % for dx = 500, moist
filename25 = 'entm_thil_dx500_moist_jul03.mat';
filename35 = 'detm_thil_dx500_moist_jul03.mat';
end

%----------------------------------------------------------------------------

% get mass flux, frac. entrainment, and detrainment matrcies
msfx1 = load(fullfile(fname,filename01));  % mass flux array  for dx=125
entrain1 = load(fullfile(fname,filename21));
detrain1 = load(fullfile(fname,filename31));

msfx = load(fullfile(fname,filename0));  % mass flux array  for dx=250
entrain = load(fullfile(fname,filename2));
detrain = load(fullfile(fname,filename3));

msfx5 = load(fullfile(fname,filename05));  % mass flux array for dx=500
entrain5 = load(fullfile(fname,filename25));
detrain5 = load(fullfile(fname,filename35));

%----------------------------------------------------------------------------
for seq = 6:6

ti = 0+(seq-1)*1; tf = 1+(seq-1)*1;

fieldname = sprintf('f%d%d',ti,tf);

mf1 = msfx1.massfx.(fieldname).mf;   % get mass flux matrix for dx=125
entm1 = entrain1.entm_thil.(fieldname).ent;  % get frac. entrainment rate matrix
detm1 = detrain1.detm_thil.(fieldname);  % get frac. detrainment rate matrix
z = msfx1.massfx.(fieldname).z;  % get z in km

mf = msfx.massfx.(fieldname).mf;   % get mass flux matrix for dx=250
entm = entrain.entm_thil.(fieldname).ent;  % get frac. entrainment rate matrix
detm = detrain.detm_thil.(fieldname);  % get frac. detrainment rate matrix
z = msfx.massfx.(fieldname).z;  % get z in km

mf5 = msfx5.massfx.(fieldname).mf;   % get mass flux matrix for dx=500
entm5 = entrain5.entm_thil.(fieldname).ent;  % get frac. entrainment rate matrix
detm5 = detrain5.detm_thil.(fieldname);  % get frac. detrainment rate matrix


%----------------------------------------------------------------------------
% Plot of z vs mass flux rate

% f = figure('units','normalized','outerposition',[0 0 1 1])
f = figure

% sp = subplot(1,3,1);
% plot(mf1,z,'-.','LineWidth',1.25,'color',[0 0.5 0]); hold on   % dx=125
% plot(mf,z,'-b','LineWidth',1.25); hold on   % dx=250
% plot(mf5,z,'--r','LineWidth',1.25); hold off   % dx=500
% 
% grid on;
% box on;
% % xlabel('mass flux \big[kgs\textsuperscript{-1}]','Interpreter','latex');
% % ylabel('Height \big[km\big]','Interpreter','latex');
% xlabel('mass flux \big[kgs\textsuperscript{-1}]','Interpreter','latex','FontSize',16);
% ylabel('Height \big[km\big]','Interpreter','latex','FontSize',16);
% 
% lgd = legend('\(\Delta=125\) m','\(\Delta=250\) m','\(\Delta=500\) m');
% set(lgd,'Interpreter','latex','Location','northeast');
% lgd.FontSize = 12;
% 
% xmin = 0;
% xmax = 3e8; 
% xtick = 1e8; 
% 
% ymin = 0;
% ymax = 10;  
% ytick = 2;
% 
% ax = gca;
% ax.XLim = [xmin,xmax];
% ax.XTick = [xmin:xtick:xmax];
% ax.XAxis.Exponent = 8;
% 
% ax.YLim = [ymin,ymax];
% ax.YTick = [ymin:ytick:ymax];
% 
% ax.FontSize = 16;
% 
% % Make y-axis:x-axis = 1.5
% pbaspect([(xmax-xmin),1.5*(xmax-xmin),1]); % multiple y-axis by the factor
% 
% tih = ti; timin = 0;
% tfh = tf; tfmin = 0;
% 
% timehi = sprintf('%02d',tih);
% timemini = sprintf('%02d',timin);
% timehf = sprintf('%02d',tfh);
% timeminf = sprintf('%02d',tfmin);
% 
% % title([num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf) ':' ...
% %     num2str(timeminf)],'Interpreter','latex');
% title([num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf) ':' ...
%     num2str(timeminf)],'Interpreter','latex','FontSize',16);
% 
% % place labels on plot
% text(0.95,1.05,'(a)','Units','normalized','Interpreter', 'latex','FontSize',20,'FontWeight','bold');
% 
% if cn == 1
%     text(0.75,0.045,'CTRL','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end
% if cn == 2
%     text(0.75,0.045,'DRY','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end
% if cn == 3
%     text(0.75,0.045,'MOIST','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end
% 
% % END OF z vs mass flux plot
% %----------------------------------------------------------------------------
%----------------------------------------------------------------------------
% % Plot of z vs frac. ENTRAINMENT rate

tol = 0.01;
entm1 = smoothing(entm1,tol);
entm = smoothing(entm,tol);
entm5 = smoothing(entm5,tol);

sp = subplot(1,2,1);
plot(entm1,z,'-.','LineWidth',1.25,'color',[0 0.5 0]); hold on   % dx=125
plot(entm,z,'-b','LineWidth',1.25); hold on   % dx=250
% plot(entmh,z,'-.','LineWidth',1.25,'color',[0.4 0.4 0.4]); hold on   % dx=250, half x-domain
plot(entm5,z,'--r','LineWidth',1.25); hold off   %dx=500
% plot(entm5h,z,':','LineWidth',1.25,'color',[0.1 0.1 0.1]); hold off   %dx=500, half x-domain

grid on;
box on;

% xlabel('\(\overline{\epsilon}\) \big[m\textsuperscript{-1}]','Interpreter','latex','FontSize',16);
xlabel('\(\overline{\epsilon}\) \big[m\(^{-1}\)]','Interpreter','latex','FontSize',16);
ylabel('Height \big[km\big]','Interpreter','latex','FontSize',16);

xmin = -6e-3;
xmax = 6e-3;
xtick = 3e-3;  % every 90 mins

ymin = 0;
ymax = 10;  
ytick = 2;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];
ax.XAxis.Exponent = -3;

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

ax.FontSize = 11;
ax.TickLabelInterpreter = 'latex';

% Make y-axis:x-axis = 1.5
pbaspect([(xmax-xmin),1.5*(xmax-xmin),1]); % multiple y-axis by the factor

tih = ti; timin = 0;
tfh = tf; tfmin = 0;

timehi = sprintf('%02d',tih);
timemini = sprintf('%02d',timin);
timehf = sprintf('%02d',tfh);
timeminf = sprintf('%02d',tfmin);

% title([num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf) ':' ...
%     num2str(timeminf)],'Interpreter','latex');
% title([num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf) ':' ...
%     num2str(timeminf)],'Interpreter','latex','FontSize',16);

% place labels on plot
text(0.88,1.05,'(a)','Units','normalized','Interpreter', 'latex','FontSize',15,'FontWeight','bold');

if cn == 1
    text(0.75,0.045,'CTRL','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
end
if cn == 2
    text(0.75,0.045,'DRY','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
end
if cn == 3
    text(0.70,0.045,'MOIST','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
end

ax1 = ax;

% 
% % END OF z vs frac. entrainment rate Plot
% %----------------------------------------------------------------------------
% 
% %----------------------------------------------------------------------------
% % Plot of z vs frac. DETRAINMENT rate

tol = 0.01*1.5;
detm1 = smoothing(detm1,tol);
detm = smoothing(detm,tol);
detm5 = smoothing(detm5,tol);
% 
% for k = 1:length(detm)
%     if detm(k) > 0.1
%         detm(k) = 0.5*(detm(k-1) + detm(k+1));
%     end
% end


% 
% tol = 0.01;
% for j = 2:length(detm5)-1
%     while abs(detm5(j)-detm5(j-1)) > tol
%         s = find(abs(detm5 - detm5(j-1)) < tol);
%         detm5(j) = 0.5*(detm5(j-1) + detm5(s(find(s>j,1))));
%     end
% end
% 
sp = subplot(1,2,2);
plot(detm1,z,'-.','LineWidth',1.25,'color',[0 0.5 0]); hold on    % dx=125
plot(detm,z,'-b','LineWidth',1.25); hold on    % dx=250
% plot(detmh,z,'-.','LineWidth',1.25,'color',[0.4 0.4 0.4]); hold on    % dx=250, half x-domain
plot(detm5,z,'--r','LineWidth',1.25); hold off    % dx=500
% plot(detm5h,z,':','LineWidth',1.25,'color',[0.1 0.1 0.1]); hold off    % dx=500, half x-domain

grid on;
box on;

xlabel('\(\overline{\delta}\) \big[m\(^{-1}\)]','Interpreter','latex','FontSize',16);
% ylabel('Height \big[km\big]','Interpreter','latex','FontSize',16);

% % lgd = legend('\(\Delta=250\) m','\(\Delta=500\) m');
% lgd = legend('\(\Delta=250\) m \(L_x=120\) km','\(\Delta=250\) m, \(L_x=60\) km',...
%     '\(\Delta=500\) m, \(L_x=120\) km','\(\Delta=500\) m, \(L_x=60\) km');
% % lgd = legend('dx250','sim75hr');
% set(lgd,'Interpreter','latex','Location','northeast');
% lgd.FontSize = 12;

xmin = -6e-3;
xmax = 6e-3;
xtick = 3e-3;  % every 90 mins

ymin = 0;
ymax = 10;  
ytick = 2;

ax = gca;
ax.XLim = [xmin,xmax];
ax.XTick = [xmin:xtick:xmax];
ax.XAxis.Exponent = -3;

ax.YLim = [ymin,ymax];
ax.YTick = [ymin:ytick:ymax];

ax.FontSize = 11;
ax.TickLabelInterpreter = 'latex';

% Make y-axis:x-axis = 1.5
pbaspect([(xmax-xmin),1.5*(xmax-xmin),1]); % multiple y-axis by the factor

tih = ti; timin = 0;
tfh = tf; tfmin = 0;

timehi = sprintf('%02d',tih);
timemini = sprintf('%02d',timin);
timehf = sprintf('%02d',tfh);
timeminf = sprintf('%02d',tfmin);

% title([num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf) ':' ...
%     num2str(timeminf)],'Interpreter','latex');
% title([num2str(timehi) ':' num2str(timemini) ' to ' num2str(timehf) ':' ...
%     num2str(timeminf)],'Interpreter','latex','FontSize',16);

% place labels on plot
text(0.88,1.05,'(b)','Units','normalized','Interpreter', 'latex','FontSize',15,'FontWeight','bold');

% if cn == 1
%     text(0.75,0.045,'CTRL','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end
% if cn == 2
%     text(0.75,0.045,'DRY','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end
% if cn == 3
%     text(0.75,0.045,'MOIST','Units','normalized','Interpreter', 'latex','FontSize',12,'FontWeight','bold');
% end

ax2 = ax;

% set ax2 left position closer to ax1
ax2.Position(1) = ax1.Position(1) + 1.15*ax1.Position(3);


% % END OF z vs DETRAINMENT plot
% %----------------------------------------------------------------------------
% %----------------------------------------------------------------------------
% 

% SAVE FIGURE

fnameg = 'C:\Users\SiuLung\Downloads';

if cn == 1
filenameg = sprintf('mfende_thil_%d%d_ctrl',ti,tf);
end

if cn == 2
filenameg = sprintf('mfende_thil_%d%d_dry',ti,tf);
end

if cn == 3
filenameg = sprintf('mfende_thil_%d%d_moist',ti,tf);
end

% saveas(sp,fullfile(fnameg,filenameg),'epsc');

end   % end seq = 6:6
