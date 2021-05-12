%% Setup
wb_dir = 'XX/workbench/bin_rh_linux64'; % set HCP Workbench directory
outdir = ''; % set directory for saving out images here

% Set parameters
num_nodes = 91282;
tr = .72; % sampling interval (s)
Fs = 1/tr; % sampling rate (Hz)
num_frames = 1191;

% HCP subjects
root_dir = 'XX/1200subject'; % HCP data directory
tasks = {'rfMRI_REST1_LR','rfMRI_REST1_RL','rfMRI_REST2_LR','rfMRI_REST2_RL'};
subjects = importdata('subjects.mat'); % subject list from Chen et al. 2020 Neuroimage

% Initialize group matrices for phase-locking values, mean network signals
% for each brain structure (cortex, striatum, thalamus, and cerebellum, and
% physiological time series
%
% (Network signals not required for computing PLVs; to save memory,
% recommended to run script separately to obtain PLVs and to obtain average
% network signals (which are used in later scripts)

plvs = single(nan(num_nodes,numel(subjects),numel(tasks)));
ctx_signals = single(nan(num_frames,8,numel(subjects),numel(tasks)));
str_signals = single(nan(num_frames,8,numel(subjects),numel(tasks)));
thal_signals = single(nan(num_frames,8,numel(subjects),numel(tasks)));
cbm_signals = single(nan(num_frames,8,numel(subjects),numel(tasks)));

physios = single(nan(num_frames,numel(subjects),numel(tasks)));


for s = 1:numel(subjects)
    tic
            
    % load BOLD time series
    subj = num2str(subjects{s});

    disp(['Processing subject ' num2str(s) ' out of ' num2str(numel(subjects))]);
    
    for t = 1:numel(tasks)

        disp(['Processing task ' num2str(t) ' out of ' num2str(numel(tasks))]);
        
        task = tasks{t};
        
        data_dir = [root_dir '/' subj '/MNINonLinear/Results/' task];
        
        try
            cifti = ft_read_cifti_mod([data_dir '/' task '_Atlas_MSMAll_hp2000_clean.dtseries.nii']);
        catch error
            disp([subj ' ' task ' missing:']);
            disp(error)
            continue
        end
        
        BOLD = cifti.data';
         
        % Physio
        try
            physio = importdata([data_dir '/' task '_Physio_log.txt']);
        catch error
            disp([subj ' ' task 'physio missing:']);
            disp(error)
            continue
        end

        % Get 6 sec windows
        if length(physio)/400<860
            continue
        end
        time_vec_bold = tr*(1:size(BOLD,1))';
        time_vec_phys = (0:length(physio)-1)'/400;
        physio_ds = zeros(size(time_vec_bold));
        
        physio(:,3) = zscore(physio(:,3));
        for i = 5:length(physio_ds)-4
            
            % For RV
            [~,phys_start] = min(abs(time_vec_phys-(time_vec_bold(i)-3)));
            [~,phys_end] = min(abs(time_vec_phys-(time_vec_bold(i)+3)));
            physio_ds(i) = std(physio(phys_start:phys_end,2));
           
            % For HRV
            %[pks,locs] = findpeaks(physio(phys_start:phys_end,3),'minpeakdistance',round(400/(180/60)));%,'minpeakwidth',400/(1/(200/60))); % max heart rate = 180 bpm; at 400 Hz, minimum of 100 samples apart
            %locs = locs(pks>prctile(physio(phys_start:phys_end,3),60));
            %Avg1(i) = mean(diff(locs))/400;
        end
        
        physio_ds = physio_ds(5:end-5);
        BOLD = BOLD(5:end-5,:);
        physio_ds = [0;diff(physio_ds)];
        

        % Get mean network signals
        if size(BOLD,1)~=1191, continue, end
        
        for n = 1:7
            ctx_signals(:,n,s,t) = nanmean(BOLD(:,mask_ctx_HCP & networks_HCP==n),2);
            str_signals(:,n,s,t) = nanmean(BOLD(:,mask_str_HCP & networks_HCP==n),2);
            thal_signals(:,n,s,t) = nanmean(BOLD(:,mask_thal_HCP & networks_HCP==n),2);
            cbm_signals(:,n,s,t) = nanmean(BOLD(:,mask_cbm_HCP & networks_HCP==n),2);
        end
            
        ctx_signals(:,8,s,t) = nanmean(BOLD(:,mask_ctx_HCP),2);
        str_signals(:,8,s,t) = nanmean(BOLD(:,mask_str_HCP),2);
        thal_signals(:,8,s,t) = nanmean(BOLD(:,mask_thal_HCP),2);
        cbm_signals(:,8,s,t) = nanmean(BOLD(:,mask_cbm_HCP),2);
        physios(:,s,t) = physio_ds;

        % Filter
        disp('Filtering')
        hp_thresh = .01; % lower bound
        lp_thresh = .05; % higher bound
        [b,a] = butter(2,[hp_thresh,lp_thresh]/(Fs/2));
        physio_ds = single(filtfilt(b,a,double(physio_ds)));
        BOLD = single(filtfilt(b,a,double(BOLD)));

        % Phase-locking values
        h1 = hilbert(physio_ds);
        h2 = hilbert(BOLD);
        plvs(:,s,t) = nanmean(exp(1i*(bsxfun(@minus,unwrap(angle(h1)),unwrap(angle(h2))))));
                
        toc
    end
end

% save('network_sigs_ctx_all.mat','ctx_signals','-v7.3')
% save('network_sigs_str_all.mat','str_signals','-v7.3')
% save('network_sigs_thal_all.mat','thal_signals','-v7.3')
% save('network_sigs_cbm_all.mat','cmb_signals','-v7.3')

%save('RV_PLVs.mat','plvs','-v7.3');