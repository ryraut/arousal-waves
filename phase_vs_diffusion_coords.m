root_dir = ''; % dir containing HCP_physio_mapping.m outputs

%% Diffusion embed

% Phase map
phases = ft_read_cifti_mod([root_dir '/HCP_RV_phasemap.dtseries.nii']);
phases = phases.data;

% FC matrix
dconn = ft_read_cifti_mod([root_dir '/HCP_S1200_1003_rfMRI_MSMAll_groupPCA_d4500ROW_zcorr.dconn.nii']);
dconn = dconn.data;

%
alpha = .5;
mask = mask_str_HCP;  % switch between str, thal, cbm
fc_mat = dconn(mask,mask_ctx_HCP);
fc_mat = bsxfun(@rdivide,fc_mat,sqrt(sum(fc_mat.^2,2)));
fc_mat = fc_mat*fc_mat';

fc_mat = 1-acos(fc_mat)/pi;
L = fc_mat;
D = sum(L,2).^-alpha;
Lalpha = L.*(D*D'); % asymmetric kernel
Dalpha = sum(Lalpha,2);
M = bsxfun(@rdivide,Lalpha,Dalpha); % graph laplacian

tic
[maps,S,~] = svd(M);
toc

%% Scatter plots
yeo_cmap = [120,18,134;70,130,180;0,118,14;196,58,250;220,248,164;230,148,34;205,62,78]/255;

figure;hold
axis off
set(gcf,'Position',[300 198 600 300])
[~,tmp1] = sort(phases(mask));
[~,tmp1] = sort(tmp1);
[~,tmp2] = sort(maps(:,2)); % principal coordinate (first nonconstant eigenvec)
[~,tmp2] = sort(tmp2);
for n=[1:4,6:7]
    scatter(tmp1(networks_HCP(mask)==n),tmp2(networks_HCP(mask)==n),20,yeo_cmap(n,:),'.')
end