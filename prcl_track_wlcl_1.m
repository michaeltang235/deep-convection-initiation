close all
clear all

%----------------------------------------------------------------------------
% REMARKS: 
% This code generates the following quantities in the chosen region each hour
% (i) mean LCL, (ii) ids of parcels found within box below the LCL at the center of domain 
% (iii) x-position of parcels found at each time step (every min.)
% (iv) z-posiiton of parcels found at each time step
% (v) ids of parcels found in region interested would then crossed the lcl
% (vi) ids of parcels which crossed the lcl and then become cloud
% (vii) ids of parcels which became cloud and then become core
% (viii) ids of parcels which became core and then acend to 6 km
% (ix) w of parcels right before crossing the lcl for the last time
%----------------------------------------------------------------------------
tic

pdata = { };

timeseq = [1:12];   % 1:12 for all hours

for seq = timeseq(1):timeseq(end)  %((ti,tf) = (0,1)) 1:12 for all hours
    
% enter range of time (hour) interested:
% ti = 0; tf = 12;
ti = 0+(seq-1)*1; tf = 1+(seq-1)*1;

dt = 4; % time step for reading cm1 output files (s=cmt seq.)
dcmt = 90; % time step between each cm1out file written (s)

% set range of x intersted:
xmin = -0.5; % initial x-pos. in km (for finding LCL)
xmax = 0.5; % final x-pos. in km

% directory to read output files (cm1 files) from:
fname = '/gs/project/qrw-161-ad/mtang/prcl_dx250_aug02';
% fname = '/storage/mtang/prcl_dx250_aug02';

% % directory for storing msfx matrix and name of matrix:
fname1 = '/gs/project/qrw-161-ad/mtang/prcl_dx250_aug02/matrices';
% fname1 = '/storage/mtang/prcl_dx250_aug02/matrices';
filename1 = sprintf('prcl_wlcl_dx250_aug02_1.mat');   % average wind vectors matrix

% enter 1 for saving output matrices and 0 for not saving outputs
saveop = 1;

%%---------------------------------------------------------------------------
%%---------------------------------------------------------------------------

% field name for structure 
fieldname = sprintf('f%d%d',ti,tf);

cmt = ti*3600/dcmt+1:dt:tf*3600/dcmt+1; %time for cm1out_00000t.nc file 

% INPUT CONSTANTS FROM constants.F
p00 = 1.0e5;
Rd = 287.04;
Cp = 1005.7;
Rv = 461.5;
epsi = Rd/Rv;
Lv = 2501000.0;

for s = cmt(1):dt:cmt(end)

filename = sprintf('cm1out_%06d.nc',s);
cmfile = fullfile(fname,filename);

% read data from output files
xf = ncread(cmfile,'xf');   % xf in km
xh = ncread(cmfile,'xh');   % xh in km
% yf = ncread(cmfile,'yf');   % yf in km
z = ncread(cmfile,'z');   % z in km

% u = ncread(cmfile,'u');      % u in m/s
% w = ncread(cmfile,'w');      % w in m/s

qv = ncread(cmfile,'qv');    % qv, water vapor in kg/kg
% qc = ncread(cmfile,'qc');    % qc, cloud water in kg/kg
% qr = ncread(cmfile,'qr');    % qr, rain 
% qs = ncread(cmfile,'qs');    % qs, snow
% qi = ncread(cmfile,'qi');    % qi, cloud ice in kg/kg
% qg = ncread(cmfile,'qg');    % qg, graupel 
% 
% qt = qv + qc + qr + qs + qi + qg; % total mixing ratio

prs = ncread(cmfile,'prs'); % pressure in Pa
th = ncread(cmfile,'th'); % potential temp. in K

% calculate temperature
exn = (prs./p00).^(Rd/Cp);       %exner function
temp = exn.*th;

% get half-level size
ni=size(qv,1);
nj=size(qv,2);
nk=size(qv,3);

% get dxh
dxh = xh(2:end) - xh(1:end-1); % in km

% get indices of xhi and xhf on xhf
xhii = max(find(abs(xh-xmin) < dxh(1)/1.5));
xhfi = min(find(abs(xh-xmax) < dxh(1)/1.5));

%-----------------------------------------------------------
% FIND LCL AT X-LOCATION INTERESTED

% convert z from km to m for parcel_free.m
zpf = 1000.*z;  % in m

% get mean (y-avg.) sounding at each x-location
% calculate y-avg. pressure, temp, and qv along x-axis
prsm = mean(prs,2);
tempm = mean(temp,2);
qvm = mean(qv,2);

% reshape prsm, tempp, and qvm into 2-d matrices as inputs
% for parcel_free.m
p_in = reshape(prsm,[ni nk]);
t_in = reshape(tempm,[ni nk]);
qv_in = reshape(qvm,[ni nk]);

% initialize LCL matrix at the beginning of time
if s == cmt(1)
    lclx = zeros(1,xhfi-xhii+1); % in m
end

% get LCL at each point on x-axis using parcel_free.m
for i = xhii:xhfi
    [~,~,lcl,~,~,~,~]=parcel_free(zpf,p_in(i,:),t_in(i,:),qv_in(i,:));
    lclx(i-xhii+1) = lclx(i-xhii+1) + lcl;
end

sprintf('cm1out_000%d.nc read, ti = %d, tf = %d \n max. lclx = %d m',s,ti,tf,max(lclx))

% get average LCL at the end of time
if s == cmt(end)
lclxm = zeros(1,xhfi-xhii+1); % in m
lclxm = lclx./(length(cmt));

sprintf('cm1out_000%d.nc read, ti = %d, tf = %d \n max. mean lclx at x intersted = %d m',s,ti,tf,max(lclxm))
end

% initialize mean lcl marix at the beginnig of time sequence
if seq == timeseq(1)
    lclxmm = zeros(length(timeseq),1);
end

% get mean lcl for each hour at the end of each cmt
if s == cmt(end)
lclxmm(seq-timeseq(1)+1) = mean(lclxm);   % in m
end

% % find index of lcl at each x in z matrix
% zid = zeros(1,length(lclx));
% for i = 1:length(lclx)
%     zid(i) = find(abs(z - lclx(i)/1000) < 0.09);   % convert lclx from m to km
% end

% sprintf('cm1out_000%d.nc read, ti = %d, tf = %d \n lclx at x intersted found',s,ti,tf)

% END FINDING LCL 
%-----------------------------------------------------------

end   % end s = cmt(1):dt:cmt(end)

%------------------------------------------------------------
% FOR pdata.nc files

% get filename
filenamep = sprintf('cm1out_pdata.nc');  % filename for parcel data
cmfilep = fullfile(fname,filenamep);

% read data from cm1out_pdata.nc
pid = ncread(cmfilep,'xh');   % xh defined as parcels id 
time = ncread(cmfilep,'time');   % time array of parcels in (s)
x = ncread(cmfilep,'x');   % x-position of parcels in (m)
y = ncread(cmfilep,'y');   % y-position of parcels in (m)
z = ncread(cmfilep,'z');   % z-position of parcels in (m)

w = ncread(cmfilep,'w');   % w velocity of parcels in (m/s)

% get mixing ratios from cm1out_pdata.nc
qv = ncread(cmfilep,'qv');    % qv, water vapor in kg/kg
qc = ncread(cmfilep,'qc');    % qc, cloud water in kg/kg
qr = ncread(cmfilep,'qr');    % qr, rain 
qs = ncread(cmfilep,'qs');    % qs, snow
qi = ncread(cmfilep,'qi');    % qi, cloud ice in kg/kg
qg = ncread(cmfilep,'qg');    % qg, graupel 

% get buoyancy from cm1out_pdata.nc
b = ncread(cmfilep,'b');   % buoyancy in m/s^2

% get qc+qi
qcqi = qc + qi;

qt = qv + qc + qr + qs + qi + qg; % total mixing ratio

% get length of pid matrix (number of parcels used)
numpr = length(pid);

% get indices of ti and tf at each cmt
tii = 1 + (seq-1)*60; 
tfi = 61 + (seq-1)*60; 

% initialize prid matrix which stores ids of parcels that are located
% within the region required (-0.5<x<0.5, -90<y<90, 0<z<lcl km)
% each hour
prid = zeros(numpr,tfi-tii+1);

% get mean lcl each hour
lclp = lclxmm(seq-timeseq(1)+1);

% region intereted: -500<x<500 (m), whole y-domain, 0<z<lcl
% find parcels in the box which have w>0
% find id of parcels that are within the region specified each hour
% parcels found are assigned a value of 1
for i = 1:numpr
    for j = tii:tfi
        if x(i,j) > -500 && x(i,j) < 500 && z(i,j) > 0 && z(i,j) < lclp && w(i,j) > 0
            prid(i,j) = 1;
        end
    end
end

% initialize pridr matrix at the beginning of time sequence
if seq == timeseq(1)
pridr = zeros(numpr,1);   % parcel id, reduced
end

% parcels found within region will have their id recorded
% record parcel's id if it is found within the region 
% update pridr if the id has not been recorded
for i = 1:numpr
    if nnz(prid(i,:)) > 0 && pridr(i) == 0  % nnz=number of nonzero entries
        pridr(i) = pid(i);
    end
end

sprintf('cm1out_pdata.nc read, ti = %d, tf = %d \n parcels located within box \n below the lcl at the center of domain found',time(tii)/3600,time(tfi)/3600) 



end  % end seq = 6:6

% END TIME SEQUENCE
%---------------------------------------------------------------------------

% parcels found within region will have their id recorded
% those not found are given value of zero
% find id of parcels with non-zero value assigned
nzprid = nonzeros(pridr);   % select non-zero entries

% find length of nzprid (number of parcels found)
lnzp = length(nzprid);

% get length of time domain
lt = size(x,2);

% get x and z positions of parcels found within the region at each time
% step
prclx = zeros(lnzp,lt);
prclz = zeros(lnzp,lt);
for i = 1:lnzp
    prclx(i,:) = x(nzprid(i),:);
    prclz(i,:) = z(nzprid(i),:);
end


%-------------------------------------------------
%-------------------------------------------------

% FIND HOW MANY PARCELS WHICH ENTERED THE REGION WOULD CROSSED THE LCL
% LATER (nlcl, flcl, and wlcl)

% initialize itlcl matrix which stores time index of parcels crossing the
% lcl for the last time
itlcl = zeros(lnzp,1);   % index of LAST time z>lcl

% initialize tlcl matrix which stores time when parcels crossing the lcl
% for the last time
tlcl = zeros(lnzp,1);

% initialize pridlcl matrix which stores id of parcels found in region interested 
% crossed the lcl later
pridlcl = zeros(lnzp,1);

% initialize itlclb matrix which stores time index of parcels right before 
% crossing the lcl for the last time
itlclb = zeros(lnzp,1);   % index of LAST time z<lcl

% initialize wlcl matrix which stores w of parcels right before crossing
% the lcl for the last time
wlcl = zeros(lnzp,1);   % in m/s

% with lcl found for each hour
% get index of time (itlcl) when parcels cross the lcl for the LAST time
% record time index, time, and parcel id if they haven't been done before
% % (avoid double counting)
% get index of time (itlclb) right before parcels crossing the lcl for the
% LAST time every hour
%
for i = 1:lnzp
    for j = timeseq(1):timeseq(end)
        if find(z(nzprid(i),1+(j-1)*60:61+(j-1)*60) > lclxmm(j-timeseq(1)+1)) ~= 0    % if z-pos. exceeds lcl of that hour
%             if itlcl(i) == 0   % if it hasn't been recorded yet
                itlcl(i) = max(find(z(nzprid(i),1+(j-1)*60:61+(j-1)*60) > lclxmm(j-timeseq(1)+1))); % index of time (in that hour) for last time crossing the lcl every hour
                tlcl(i) = time(1+(j-1)*60-1+itlcl(i)); % get time from time index + index of j th hour
                itlclb(i) = max(find(z(nzprid(i),1:tlcl(i)/60+1) < lclxmm(j-timeseq(1)+1)));   % index of time right before crossing lcl last time every hour
%                 if z(nzprid(i),tlcl(i)/60+1) >= z(nzprid(i),tlcl(i)/60) 
            if time(itlclb(i)) < tlcl(i)
                    pridlcl(i) = nzprid(i);
                    wlcl(i) = w(nzprid(i),itlclb(i));   % update wlcl every hour
%                 else
%                     itlcl(i) = 0;
%                     tlcl(i) = 0;
            end
%             end
%                 itlclb(i) = max(find(z(nzprid(i),1:tlcl(i)/60+1) < lclxmm(j-timeseq(1)+1)));   % index of time right before crossing lcl last time every hour
%                 wlcl(i) = w(nzprid(i),itlclb(i));   % update wlcl every hour
        end
    end
end

% get id of parcels in pridlcl with non-zero value assigned
nzpridlcl = nonzeros(pridlcl);   % select non-zero entries

% get w of parcels in wlcl with non-zero value assigned
wlclr = nonzeros(wlcl);   % select non-zero entries

% get N = number of parcels enter the region interested
N = lnzp;

% get Nlcl = number of parcels found in the region interested crossed the
% LCL later
Nlcl = length(nzpridlcl);

% get flcl which is fraction of Nlcl/N
flcl = Nlcl/N;

% get x and z positions of parcels which were found within the region and would cross the lcl later 
% at each time step
prclxlcl = zeros(Nlcl,lt);
prclzlcl = zeros(Nlcl,lt);
for i = 1:Nlcl
    prclxlcl(i,:) = x(nzpridlcl(i),:);
    prclzlcl(i,:) = z(nzpridlcl(i),:);
end

sprintf('flcl = %d',flcl) 


% END FINDING HOW MANY PARCELS IN THE REGION INTERESTED CROSSED THE LCL (Nlcl, flcl, wlcl)
%-------------------------------------------------------------
%-------------------------------------------------------------

% FIND HOW MANY OF PARCELS CORSSED THE LCL WOULD BECOME CLOUD 
% threhold qth: qc+qi>0.1 g/kg

% set threshold qth
qth = 0.1;  % threshold qc+qi in g/kg  
qth = qth/1000;  % qth in kg/kg

% initialize prid_qcqi matrix which stores id of parcels which have crossed the LCL and have qc+qi>0.1
% g/kg at anytime, those not are assigned value of zero
pridqcqi = zeros(Nlcl,1);
% for i = 1:Nlcl
%     if find(qcqi(nzpridlcl(i),:) > qth) ~= 0
%         pridqcqi(i) = nzpridlcl(i);
%     end
% end

% get number of parcels from N_lcl would become cloud
% parcels have to have qc+qi>0.1 g/kg and w>0 at the same time
for i = 1:Nlcl
    for j = 1:lt
        if qcqi(nzpridlcl(i),j) > qth && w(nzpridlcl(i),j) > 0
            pridqcqi(i) = nzpridlcl(i);
        end
    end
end

% find id of parcels in prid_qcqi with non-zero value assigned
nzpridqcqi = nonzeros(pridqcqi);   % select non-zero entries

% find number of parcels with qc+qi>0.1 g/kg among those already found
% within the region and have already passed the lcl
% find length of nzpridqcqi 
Ncloud = length(nzpridqcqi); 

% get fcloud which is fraction of Ncloud/Nlcl
fcloud = Ncloud/Nlcl;

sprintf('fcloud = %d',fcloud) 

% 
% % END FIND HOW MANY PARCELS FROM N_LCL BECOME CLOUD (HAVE qc+qi>0.1 g/kg) 
% %-------------------------------------------------
% %-------------------------------------------------

% FIND HOW MANY PARCELS WHICH CROSSED THE LCL AND BECOME CLOUD LATER WOULD
% THEN BECOME CORE

% initialize prid_wbpqcqi matrix which stores id of parcels from Ncloud that would become core, 
% those not are assigned value of zero
pridwbpqcqi = zeros(Ncloud,1);

% get number of parcels from N_cloud that would become core
% parcels have to have w>0, b>0, and qc+qi>0.1 g/kg at the same time
% only record such parcels once
for i = 1:Ncloud
    for j = 1:lt
        if w(nzpridqcqi(i),j) > 0 && b(nzpridqcqi(i),j) > 0 && qcqi(nzpridqcqi(i),j) > qth && pridwbpqcqi(i) == 0
            pridwbpqcqi(i) = nzpridqcqi(i);
        end
    end
end

% find id of parcels in prid_wbpqcqi with non-zero value assigned
nzpridwbpqcqi = nonzeros(pridwbpqcqi);   % select non-zero entries

% find number of parcels having w>0, b>0, and qc+qi>0.1 g/kg among those that have already 
% been found within the region and have already been found to have passed the lcl and have become cloud
% find length of nzpridwbpqcqi 
Ncore = length(nzpridwbpqcqi); 

% get fcloud which is fraction of Ncore/Ncloud
fcore = Ncore/Ncloud;

sprintf('fcore = %d',fcore) 

% % END FIND HOW MANY PARCELS FROM N_CLOUD BECOME CORE (HAVE w>0, b>0, and qc+qi>0.1 g/kg) 
% %--------------------------------------------------------------
% %--------------------------------------------------------------

% FIND HOW MANY PARCELS WHICH BECAME CLOUD WOULD ASCEND TO 6km 

% initialize prid_coreha matrix which stores id of parcels from Ncore that
% would ascend to 6 km 
pridcoreha = zeros(Ncore,1);   % parcel id, core of high altitude

% get number of parcels from N_core that would acend to 6 km
% only record such parcels once
for i = 1:Ncore
    if find(z(nzpridwbpqcqi(i),:) > 6000) ~= 0
        pridcoreha(i) = nzpridwbpqcqi(i);
    end
end

% find id of parcels in prid_core_ha with non-zero value assigned
nzpridcoreha = nonzeros(pridcoreha);   % select non-zero entries

% find number of parcels acended to 6 km among those that have already 
% been found determined as core
% find length of nzpridcoreha
Ncoreha = length(nzpridcoreha); 

% get fcloud which is fraction of Ncore_ha/Ncore
fcoreha = Ncoreha/Ncore;

sprintf('fcoreha = %d',fcoreha) 

% % END FIND HOW MANY PARCELS FROM N_CORE BECOME ACEND TO 6 km
% %--------------------------------------------------------------
% %--------------------------------------------------------------
    
% 
% % store matrices to pdata array
pdata.lclxmm = lclxmm;    % mean of LCL over region interested every hour
pdata.nzprid = nzprid;   % id of parcels found within region
pdata.prclx = prclx;   % x-position (m) of parcels found within region 
pdata.prclz = prclz;   % z-position (m) of parcels found within region
pdata.prclxlcl = prclxlcl;   % x-position (m) of parcels found within region that would cross the lcl
pdata.prclzlcl = prclzlcl;   % z-position (m) of parcels found within region that would cross the lcl
pdata.nzpridlcl = nzpridlcl;   % id of parcels found in the region that would cross the lcl 
pdata.nzpridqcqi = nzpridqcqi;   % id of parcels which crossed the lcl and then become cloud
pdata.nzpridwbpqcqi = nzpridwbpqcqi;   % id of parcels which became cloud and then become core
pdata.nzpridcoreha = nzpridcoreha;   % id pf parcels which became core and then acend to 6 km
pdata.flcl = flcl;   % flcl = Nlcl/N, N = no. of parcels entered region interested
pdata.fcloud = fcloud;   % fcloud = Ncloud/Nlcl
pdata.fcore = fcore;   % fcloud = Ncore/Ncloud
pdata.fcoreha = fcoreha;   % fcloudha = Ncoreha/Ncloud, Ncoreha = cores that reach 6 km
pdata.wlclr = wlclr;   % w (in m/s) of parcels right before crossing the lcl for the last time
pdata.xf = xf; % in km
pdata.z = z; % in km
%  

% save workspace variables 
if saveop == 1
save(fullfile(fname1,filename1),'pdata');
end
    

toc

% pdata.nzpridbp = nzpridbp;   % id of parcels found within region with b>0
% pdata.nzpridqcqi = nzpridqcqi;   % id of parcels found within region with qc+qi>0.1 g/kg
% pdata.nzpridbpqcqi = nzpridbpqcqi;   % id of parcels found within region with b>0 and qc+qi>0.1 g/kg
% pdata.rbp = rbp;   % ratio of parcels in region interested with b>0.1 m/s^2
% pdata.rqcqi = rqcqi;   % ratio of parcels in region interested with qc+qi>0.1 g/kg
% pdata.rbpqcqi = rbpqcqi;   % ratio of parcels in region interested with b>0 m/s^2 and qc+qi>0.1 g/kg
    
    



