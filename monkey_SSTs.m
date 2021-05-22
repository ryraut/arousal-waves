%% Inputs
monkey = 'Chibi'; % Chibi, George, Kin2, Su
cond = 2; % 1 = EO; 2 = EC; 3 = sleep; 4 = anesth

%% Setup
root_dir = ['/data/nil-bluearc/raichle/ryan/ecog/' monkey];
load([root_dir '/' monkey 'Map.mat']);

switch cond
    case 1
        label1 = 'EO';
        label2 = 'AwakeEyesOpened';
    case 2
        label1 = 'EC';
        label2 = 'AwakeEyesClosed';
    case 3
        label1 = 'sleep';
        label2 = 'Sleeping';
    case 4
        label1 = 'anesth';
        label2 = 'Anesthetized';
end

runs = textread([root_dir '/' monkey '_' label1 '_runs.txt'],'%s');

num_nodes = 128;
Fs = 1000;

% Initialize SST variables
avg_sst = zeros(201,101);
avg_gamma = zeros(201,128);
avg_alpha = avg_gamma;
avg_delta = avg_gamma;
avg_acfs = avg_gamma;
sst_count = 0;
null_sst = zeros(201,101);
null_count = 0;

for r = 1:length(runs)
    disp(['Processing ' runs{r}])
    disp([num2str(r) ' out of ' num2str(numel(runs))])
    data_dir = [root_dir '/../' runs{r}];
    time_vec = importdata([data_dir '/ECoGTime.mat']);
    info = importdata([data_dir '/Condition.mat']);
    start = info.ConditionTime(strcmp(info.ConditionLabel,[label2 '-Start'])); % in seconds
    stop = info.ConditionTime(strcmp(info.ConditionLabel,[label2 '-End']));

    start = find(time_vec==start);
    if isempty(start)
        start = 10*Fs; % start ten seconds in
    end
    stop = find(time_vec==stop);
    if isempty(stop)
        stop = length(time_vec);
    end
    num_samples = length(start:stop);

    %% Load data
    disp('Loading data...')
    raw_data = single(zeros(num_samples,num_nodes));
    for i = 1:num_nodes
        tic
        temp = importdata([data_dir '/ECoG_ch' num2str(i) '.mat']);
        raw_data(:,i) = temp(start:stop);
        toc
    end

    good = true(num_nodes,1);
    if strcmp(monkey,'George')
        good([53,73]) = false;
    elseif strcmp(monkey,'Su')
        good(50) = false;
    end
    raw_data(:,~good) = nan;
    raw_data = bsxfun(@minus,raw_data,nanmean(raw_data,2));
    
    %% Time-frequency analysis
    
    [b1,a1] = butter(4,[49,51]/(Fs/2),'stop'); % notch
    [b2,a2] = butter(4,[99,101]/(Fs/2),'stop'); % notch
    filt_data = filtfilt(b1,a1,double(raw_data));
    filt_data = filtfilt(b2,a2,double(filt_data));
        
    clear struct
    Fs_smooth = 5;
    struct.Fs = 1000;
    struct.fpass = [1 100];
    struct.trialave = 0;
    struct.tapers = [3 5];
    [S,t,f] = mtspecgramc(filt_data,[1,1/Fs_smooth],struct);

    power = real(10*log10(S));
    power_zm = bsxfun(@minus,power,nanmean(power,1));
    global_power = nanmean(power_zm,3)';
    
%     figure;
%     surf(t,f,global_power,'EdgeColor','none');
%     axis xy;axis tight;colormap(jet);view(0,90);caxis([-5,5])
%     colorbar;c=colorbar;c.Label.String = 'Power (dB)';
%     xlabel('Seconds')
%     ylabel('Frequency (Hz)')

    %% SST detection

    % Get low-frequency power
    [~,ind_low] = min(abs(f-0));
    [~,ind_high] = min(abs(f-4));
    lf_power = squeeze(nanmean(power_zm(:,ind_low:ind_high,:),2));

    % Get mid-frequency power
    [~,ind_low] = min(abs(f-9));
    [~,ind_high] = min(abs(f-21));
    mf_power = squeeze(nanmean(power_zm(:,ind_low:ind_high,:),2));

    % Get gamma BLP
    [~,ind_low] = min(abs(f-42));
    [~,ind_high] = min(abs(f-87));
    hf_power = squeeze(nanmean(power_zm(:,ind_low:ind_high,:),2));
    
    % Find SSTs
    stdevs = nanstd(lf_power);
    means = nanmean(lf_power);
    sst_inds = sum(bsxfun(@gt,lf_power,(means+stdevs)),2)>(.4*sum(good));
   
    [b,a] = butter(2,[.01,.05]/(Fs_smooth/2)); % infra-slow envelope of gamma blp
    hf_power = filtfilt(b,a,hf_power);
    mf_power = filtfilt(b,a,mf_power);
    lf_power = filtfilt(b,a,lf_power);


    last_sst = -100;
    for i = 1:length(sst_inds)
        if sst_inds(i) && i>101 && i<(length(sst_inds)-100) && (i-3*Fs_smooth)>last_sst
            avg_sst = avg_sst + global_power(:,i-100:i+100)';
            avg_gamma = avg_gamma + hf_power(i-100:i+100,:);
            avg_alpha = avg_alpha + mf_power(i-100:i+100,:);
            avg_delta = avg_delta + lf_power(i-100:i+100,:);

            sst_count = sst_count+1;
            last_sst = i;
        end
    end

    % Generate null distribution
    last_sst = -100;
    null_inds = sst_inds(randperm(length(sst_inds)));
    for i = 1:length(null_inds)
        if null_inds(i) && i>101 && i<(length(null_inds)-100) && (i-3/.2)>last_sst && null_count <= sst_count
            null_sst = cat(3,null_sst,global_power(:,i-100:i+100)');
            null_count = null_count+1;
            last_sst = i;
        end
    end
end

avg_sst = avg_sst./sst_count;
avg_gamma = avg_gamma./sst_count;
avg_alpha = avg_alpha./sst_count;
avg_delta = avg_delta./sst_count;

% normalize by null
avg_sst = avg_sst-squeeze(mean(null_sst,3));
avg_sst = avg_sst./std(null_sst,0,3);

figure;
surf((1/Fs_smooth)*(-100:100),f,avg_sst','EdgeColor','none');
axis xy;axis tight;colormap(jet);view(0,90);
colorbar;c=colorbar;%c.Label.String = 'Z';
xlabel('Seconds')
ylabel('Frequency (Hz)')
caxis([-1,1])
set(gca,'fontsize',13,'fontweight','bold')

figure;plot(nanmean(avg_gamma,2))
