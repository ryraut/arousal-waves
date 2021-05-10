% Diffusion embed
phases = ft_read_cifti_mod([root_dir '/HCP_RV_phasemap.dtseries.nii']);
phases = phases.data;

dconn = ft_read_cifti_mod([root_dir '/HCP_S1200_1003_rfMRI_MSMAll_groupPCA_d4500ROW_zcorr.dconn.nii']);
dconn = dconn.data;

%
alpha = .5;
mask1 = mask_str_HCP;
fc_mat = dconn(mask1,mask_ctx_HCP);
fc_mat = bsxfun(@rdivide,fc_mat,sqrt(sum(fc_mat.^2,2)));
fc_mat = fc_mat*fc_mat';

fc_mat = 1-acos(fc_mat)/pi;
L = fc_mat;
D = sum(L,2).^-alpha;
Lalpha = L.*(D*D');
Dalpha = sum(Lalpha,2);
M = bsxfun(@rdivide,Lalpha,Dalpha);
[maps,S,~] = svd(M);

figure;hold
axis off
set(gcf,'Position',[300 198 600 300])
[~,tmp1] = sort(phases(mask1));
[~,tmp1] = sort(tmp1);
[~,tmp2] = sort(maps(:,2));
[~,tmp2] = sort(tmp2);
for n=[1:4,6:7]
    scatter(tmp1(networks_HCP(mask1)==n),tmp2(networks_HCP(mask1)==n),20,yeo_cmap(n,:),'.')
end


% Scatter plots
figure;hold
axis off
set(gcf,'Position',[300 198 600 300])
[~,tmp1] = sort(temp1(mask_ctx_HCP));
[~,tmp1] = sort(tmp1);
[~,tmp2] = sort(-grads(mask_ctx_HCP));
[~,tmp2] = sort(tmp2);
for n=1:7
    scatter(tmp1(networks_HCP(mask_ctx_HCP)==n),tmp2(networks_HCP(mask_ctx_HCP)==n),5,yeo_cmap(n,:),'.')
end