% import files
root_dir = ''; % directory containing outputs from "main"
physio = -reshape(importdata([root_dir '/physios_all.mat']),1191,[]);
globals = importdata([root_dir '/network_sigs_all.mat']);globals=reshape(globals(:,8,:,:),1191,[]);
nanmask = sum(isnan(physio)) | sum(isnan(globals));
nanmask = nanmask | sum(physio)==0 | sum(globals)==0;
physio = detrend(physio(:,~nanmask));
globals = detrend(globals(:,~nanmask));

% Coherence
clear struct
tr = .72; % sampling interval (s)
Fs = 1/tr; % sampling rate (Hz)

% Define params for coherence
struct.Fs = Fs;
struct.fpass = [0 .12];
struct.trialave = 1;
struct.tapers = [6 6];

% Compute global signal phase spectrum (to subtract out below)
d = round(7/tr); % initial time shift to faciliate phase shift measurement in presence of long delays (probably unnecessary)
globals2 = [globals(d+1:end,:);zeros(d,size(globals,2));];
[C,phi_global,S12,S1,S2,f] = coherencyc(globals2,physio,struct);

figure;hold
set(gca,'fontsize',20,'fontweight','bold')
ax = gca;
ax.LineWidth = 2;
for n = [2,4,6,7]
    tic
    net_sigs = importdata([root_dir '/network_sigs_all.mat']);net_sigs=reshape(net_sigs(:,n,:,:),1191,[]);
    net_sigs = net_sigs(:,~nanmask);
    d = round(7/tr);
    net_sigs = [net_sigs(d+1:end,:);zeros(d,size(net_sigs,2));]; % 7 sec time shift
    [C,phi,S12,S1,S2,f] = coherencyc(net_sigs,physio,struct);
    %plot(f,C,'color',yeo_cmap(n,:),'linewidth',3)
    plot(f,unwrap(phi-phi_global),'color',yeo_cmap(n,:),'linewidth',3)
    toc
end

% create null coherence values
num_trials = 500;
net_sigs = nanmean(importdata([root_dir '/network_sigs_all.mat']),2);
net_sigs=reshape(net_sigs,1191,[]);net_sigs = net_sigs(:,~nanmask);
C_nulls = zeros(length(f),num_trials);
for n = 1:num_trials
    disp(num2str(n))
    [temp,~,~,~,~,~] = coherencyc(net_sigs,physio(:,randperm(size(physio,2))),struct); % shuffle physiological time series across subjects
    C_nulls(:,n) = nanmean(temp,2);
end
C_null = prctile(C_nulls,99,2);
plot(f,C_null,'color','k','linewidth',2,'linestyle','--')
axis([.01,.1,0,.5])

% phase
ax.YTickLabel = {'-\pi','-\pi/2','0','\pi/2','\pi'};
ax.YTick = [-pi -pi/2 0 pi/2 pi];
axis([.01,.1,-pi,pi])