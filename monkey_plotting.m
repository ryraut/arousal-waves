% Plot SST spectrogram and gamma band time series

SST1 = importdata('XX/Chibi/SST_EC.mat');
avg_sst1 = SST1.avg_sst.*std(SST1.null_sst,0,3)+squeeze(mean(SST1.null_sst,3));
sst_count1 = SST1.sst_count;
avg_sst1 = avg_sst1.*sst_count1;

SST2 = importdata('XX/George/SST_EC.mat');
avg_sst2 = SST2.avg_sst.*std(SST2.null_sst,0,3)+squeeze(mean(SST2.null_sst,3));
sst_count2 = SST2.sst_count;
avg_sst2 = avg_sst2.*sst_count2;
avg_sst = (avg_sst1+avg_sst2)./(sst_count1+sst_count2);
Fs_smooth = 5;
f = importdata('XX/George/f_vec.mat');

figure;
surf((1/Fs_smooth)*(-100:100),f,avg_sst','EdgeColor','none');
axis xy;axis tight;colormap(viridis);view(0,90);
colorbar;caxis([-.25,.25])
set(gca,'fontsize',15,'fontweight','bold')


% Plot lags
%avg_gamma = SST1.avg_gamma;
%load('XX/Chibi/ChibiMap.mat');

avg_gamma = SST2.avg_gamma;
load('XX/George/GeorgeMap.mat');

lag_lim = 7;
Fs_smooth = 5;
num_nodes = size(avg_gamma,2);
data = avg_gamma;
data1 = nanmean(data,2);
data2 = data;
data1 = data1 - nanmean(data1);
data2 = bsxfun(@minus,data2,nanmean(data2));

ccfs = zeros(1+2*lag_lim*Fs_smooth,num_nodes);
for i = 1:num_nodes
    tic
    ccfs(:,i) = xcorr(data2(:,i),data1,lag_lim*Fs_smooth,'coeff');
    toc
end

% store
[corrs,inds] = max(ccfs);
inds(isnan(corrs)) = nan;
inds = inds-(length(ccfs)-1)/2;
inds = inds-nanmedian(inds);
inds = inds/Fs_smooth;

figure
image(I);axis equal
hold on
scatter(X,Y,50,inds,'filled','linewidth',2);
axis off


% Plot time series
figure;hold; box off
set(gca,'fontsize',13,'fontweight','bold')
jets = jet(128);
[a,b] = sort(inds);
avg_gamma_sort = avg_gamma(:,b);
for n = 1:128
    p1 = plot(-20:1/Fs_smooth:20,avg_gamma_sort(:,n),'color',jets(n,:),'linewidth',1);
    p1.Color(4) = .75;
end
plot(-20:1/Fs_smooth:20,nanmean(avg_gamma,2),'k','linewidth',5)
axis off

obj = scalebar;
obj.YLen = .1;
obj.XLen = 5;
obj.YUnit = 'dB';
obj.XUnit = 's';
obj.hTextX_Pos = [1.4 -.012];
obj.Position = [-18 -.18];