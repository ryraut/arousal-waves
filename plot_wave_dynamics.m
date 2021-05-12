% Phase-amplitude plot

% import files
physios = reshape(importdata('RV_all.mat'),1191,[]);
globals = importdata('ctx_sigs_all.mat');globals=reshape(globals(:,8,:,:),1191,[]);
nanmask = sum(isnan(physios)) | sum(isnan(globals));
physios = physios(:,~nanmask);

% Filter
hp_thresh = .01;   % lower bound
lp_thresh = .05;   % upper bound
tr = .72;
Fs = 1/tr;
[b,a] = butter(2,[hp_thresh,lp_thresh]/(Fs/2));
physios = filtfilt(b,a,double(physios));
num_bins = 40;
bin_lims = linspace(-pi,pi,num_bins+1)';
bins = zeros(num_bins,size(physios,2));
smooth_width = 3;
wrap_pupils = angle(hilbert(physios));

yeo_cmap = [120,18,134;70,130,180;0,118,14;196,58,250;220,248,164;230,148,34;205,62,78]/255;
figure;hold
set(gca,'fontsize',18,'fontweight','bold')
set(gcf,'Position',[300 198 650 315])
midpoints = conv(bin_lims,[.5,.5],'valid');

for RSN=[2,3,4,6,7]
    net_sigs = importdata('cbm_signals_all.mat');net_sigs=reshape(net_sigs(:,RSN,:,:),1191,[]);
    net_sigs = detrend(net_sigs(:,~nanmask));
    
    % Filter
    hp_thresh = .01;   % lower bound
    lp_thresh = .05;   % upper bound
    [b,a] = butter(1,[hp_thresh,lp_thresh]/(Fs/2));
    net_sigs = filtfilt(b,a,double(net_sigs));
    
    disp('Computing real')    
    for b = 2:num_bins+1
        net_sigs2 = net_sigs;net_sigs2(wrap_pupils<bin_lims(b-1) | wrap_pupils>bin_lims(b))=nan;
        bins(b-1,:) = nanmean(net_sigs2,1);
    end

    % smooth bins
    bins3 = bins;
    bins2 = repmat(bins,3,1);
    for i = num_bins+1:2*num_bins
        bins3(i-num_bins,:) = nanmean(bins2(i-smooth_width:i+smooth_width,:),1);
    end

    shadedErrorBar(midpoints,nanmean(bins3,2),(1/sqrt(size(physios,2)))*nanstd(bins3,0,2),'lineProps',{'color',yeo_cmap(RSN,:)})
    plot(midpoints,nanmean(bins3,2),'color',yeo_cmap(RSN,:),'linewidth',2)

end

% generate null distribution
disp('Generating null')
ntrials = 200;
nulls = zeros(num_bins,ntrials);
net_sigs = nanmean(importdata('cbm_sigs_all.mat'),2);net_sigs=reshape(net_sigs,1191,[]);
net_sigs = detrend(net_sigs(:,~nanmask));
for n = 1:ntrials
    disp(n)
    rand_phases = wrap_pupils(:,randperm(size(wrap_pupils,2)));
    
    for b = 2:num_bins+1
        net_sigs2 = net_sigs;net_sigs2(rand_phases<bin_lims(b-1) | rand_phases>bin_lims(b))=nan;
        bins(b-1,:) = nanmean(net_sigs2,1);
    end
    
    % Smooth bins
    bins2 = bins;
    bins3 = bins2;
    bins2 = repmat(bins2,3,1);
    for i = num_bins+1:2*num_bins
        bins3(i-num_bins,:) = nanmean(bins2(i-smooth_width:i+smooth_width,:),1);
    end
     
    nulls(:,n) = nanmean(bins3,2);
end

shadedErrorBar(midpoints,mean(nulls,2),[prctile(nulls,99,2)-mean(nulls,2),mean(nulls,2)-prctile(nulls,1,2)]','lineProps','k')
plot(midpoints,mean(nulls,2),'k','linewidth',2)

amp = 14;
plot(-pi:.01:0,nanmean(bins3(:))+amp*cos(-pi:.01:0),'r','linewidth',2) % .02 for FD
plot(0:.01:pi,nanmean(bins3(:))+amp*cos(0:.01:pi),'b','linewidth',2)

axis([-pi,pi,-14,14.3])
ax = gca;
ax.XTick = [-pi -pi/2 0 pi/2 pi];
ax.XTickLabel = {'-\pi','-\pi/2','0','\pi/2','\pi'};
ax.YTickLabel = {''};
ax.YTick = 0;


