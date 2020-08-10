
function mfmf = fracmf(massfx,z,cmt)

dz = 1000.*(z(2:end) - z(1:end-1));  % dz in m

ni = size(massfx,1);
nj = size(massfx,2);

%create dmf/dz matrix
dmfdz = zeros(ni,nj);

% % central difference for massfx(2:ni-1,:)
for i = 2:ni-1
    dmfdz(i,:) = (massfx(i+1,:) - massfx(i-1,:))./(dz(i) + dz(i-1));
end

% % forward difference for massfx(1,:):
dmfdz(1,:) = (massfx(2,:) - massfx(1,:))./dz(1);

% % backward difference for massfx(ni,:):
dmfdz(ni,:) = (massfx(ni,:) - massfx(ni-1,:))./dz(end);

% calculate fractional change in dmf/dz:
fmf = zeros(ni,nj);

for i = 1:ni
    for j = 1:nj
        if massfx(i,j) == 0
            fmf(i,j) = 0;
        else
            fmf(i,j) = dmfdz(i,j)/massfx(i,j);
        end
    end
end

% mean fractional dmf/dz
mfmf = mean(fmf(:,cmt(1):cmt(end)),2);

for i = 1:ni
    if mfmf(i) < -0.08
       mfmf(i) = 0.5*(mfmf(i-1) + mfmf(i+1));
    end
    if mfmf(i) > 0.08
        mfmf(i) = 0.5*(mfmf(i-1) + mfmf(i+1));
    end
    
end
