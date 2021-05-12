%% Flow fields

% 3D extension of https://www.mathworks.com/matlabcentral/fileexchange/22756-horn-schunck-optical-flow-method
% makes use of CIFTI tools

wbdir = 'XX/workbench/bin_rh_linux64';
cifti = ft_read_cifti_mod('RV_dynamics_all_sm3_zm.dtseries.nii');
movie = cifti.data';
left_mask = cifti.brainstructure(cifti.brainstructure>0)==1;
brain_mask = mask_ctx_HCP;
system([wbdir '/wb_command -cifti-gradient RV_dynamics_all_sm3_zm.dtseries.nii COLUMN temp_gradients.dtseries.nii -left-surface ' surf_dir '/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.very_inflated.32k_fs_LR.surf.gii -right-surface ' surf_dir '/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.very_inflated.32k_fs_LR.surf.gii -surface-presmooth 20 -volume-presmooth 10 -vectors temp_vecs.dscalar.nii'])
vecs32k = ft_read_cifti_mod('test_vecs.dscalar.nii');
vecs32k = vecs32k.data(brain_mask,:);
vecs32k = reshape(vecs32k,sum(brain_mask),3,[]);
FX = squeeze(vecs32k(:,1,2:end))';
FY = squeeze(vecs32k(:,2,2:end))';
FZ = squeeze(vecs32k(:,3,2:end))';
FT = diff(movie(:,brain_mask));

% Iterations
UX = zeros(sum(brain_mask),39);
UY = UX;
UZ = UX;
ux = zeros(size(FX,2),1);
uy = ux;
uz = ux;
iters = 1000;
alpha = 5;
kernel = zeros(3,3,3);
kernel(:,:,1) = [0 0 0;0 1 0;0 0 0]/6;
kernel(:,:,2) = [0 1 0;1 0 1;0 1 0]/6;
kernel(:,:,3) = kernel(:,:,1);
for frame = 1:size(UX,2)
    tic
    fx = FX(frame,:)';
    fy = FY(frame,:)';
    fz = FZ(frame,:)';
    ft = FT(frame,:)';
    for i=1:iters
        % Compute local averages of the flow vectors
        uxAvg=convn(ux,kernel,'same');
        uyAvg=convn(uy,kernel,'same');
        uzAvg=convn(uz,kernel,'same');
        % Compute flow vectors constrained by its local average and the optical flow constraints
        ux=uxAvg - ( fx.*( (fx.*uxAvg) + (fy.*uyAvg) + (fz.*uzAvg) + ft))...
        ./ ( alpha.^2 + fx.^2 + fy.^ 2 + fz.^ 2);
        uy=uyAvg - ( fy.*( (fx.*uxAvg) + (fy.*uyAvg) + (fz.*uzAvg) + ft))...
        ./ ( alpha.^2 + fx.^2 + fy.^ 2 + fz.^ 2);
        uz=uzAvg - ( fz.*( (fx.*uxAvg) + (fy.*uyAvg) + (fz.*uzAvg) + ft))...
        ./ ( alpha.^2 + fx.^2 + fy.^ 2 + fz.^ 2);
    end
    UX(:,frame) = ux;
    UY(:,frame) = uy;
    UZ(:,frame) = uz;
    toc
end
UX1 = nanmean(UX,2);UY1 = nanmean(UY,2);UZ1 = nanmean(UZ,2);

% for ctx -- downsample to 4k and plot
modes = [UX1,UY1,UZ1];
cifti = ft_read_cifti_mod('RV_dynamics_all_sm3_zm.dtseries.nii');
cifti.hdr.dim(6)=3;cifti.data=zeros(size(cifti.data,1),3);
cifti.data(brain_mask,:)=modes;
ft_write_cifti_mod('temp.dtseries.nii',cifti);
system([wbdir '/wb_command -cifti-resample temp.dtseries.nii COLUMN ' surf_dir '/template.4k.dtseries.nii COLUMN ADAP_BARY_AREA CUBIC temp.4k.dtseries.nii -surface-largest -left-spheres ' surf_dir '/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.sphere.32k_fs_LR.surf.gii ' surf_dir '/Sphere.4k.L.surf.gii -left-area-surfs ' surf_dir '/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.midthickness.32k_fs_LR.surf.gii ' surf_dir '/Conte69.L.midthickness.4k_fs_LR.surf.gii -right-spheres ' surf_dir '/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.sphere.32k_fs_LR.surf.gii ' surf_dir '/Sphere.4k.R.surf.gii -right-area-surfs ' surf_dir ' /Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.midthickness.32k_fs_LR.surf.gii ' surf_dir '/Conte69.R.midthickness.4k_fs_LR.surf.gii']);
temp = ft_read_cifti_mod('temp.4k.dtseries.nii');
coords_4k = gifti([surf_dir '/Conte69.L.very_inflated.4k_fs_LR.surf.gii']);
mask = gifti([surf_dir '/nonmedial_wall.L.4k.func.gii']);
coords_4k = coords_4k.vertices(mask.cdata==1,:);
x=coords_4k(:,1);y=coords_4k(:,2);z=coords_4k(:,3);
cifti = ft_read_cifti_mod('RV_dynamics_all_sm3_zm.4k.dtseries.nii');
brain_mask_4k = cifti.brainstructure(cifti.brainstructure>0)==1;
UX1=temp.data(brain_mask_4k,1);UY1=temp.data(brain_mask_4k,2);UZ1=temp.data(brain_mask_4k,3);
modes = [UX1,UY1,UZ1];
coords = gifti([surf_dir '/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.very_inflated.32k_fs_LR.surf.gii']);
coords = coords.vertices;

inds = ceil(atan2(UY1,UZ1)/(pi/180)); % left lateral
%inds = ceil(atan2(-UY1,UZ1)/(pi/180)); % left medial
inds(inds<1) = inds(inds<1)+360;


% colorwheel
deg=linspace(0,2*pi,21);
u=cos(deg+180)';
v=sin(deg+180)';
inds=ceil(atan2(-u,v)/(pi/180));
inds(inds<1)=inds(inds<1)+360;
figure;hold;axis equal;cmap = colormap(hsv(360));axis off;set(gcf,'Color','k')
for i=1:360
    quiver(zeros(sum(inds==i),1),zeros(sum(inds==i),1),u(inds==i),v(inds==i),'color',cmap(i,:),'linewidth',8,'maxheadsize',.5)
end
scatter(0,0,500,'k','filled')



% compute surface normals
norms = gifti('L.very_inflated_normals.4k_fs_LR.func.gii');
mask_4k = gifti('/data/nil-bluearc/raichle/ryan/surface/nonmedial_wall.L.4k.func.gii');
norms = norms.cdata(mask_4k.cdata==1,:);

% for lateral view
inds = ceil(atan2(UY1,UZ1)/(pi/180));
inds(inds<1) = inds(inds<1)+360;
inds(norms(:,1)>0)=nan;

% for medial view
inds = ceil(atan2(-UY1,UZ1)/(pi/180));
inds(inds<1) = inds(inds<1)+360;
inds(norms(:,1)<0)=nan;
inds(x<0) = nan;

% plot
figure;hold;axis off
set(gcf,'Position',[300 198 650 420])
set(gcf,'color','k');
%scatter3(coords(:,1),coords(:,2),coords(:,3),10,'k','filled')
cmap = colormap(hsv(360));
for i = 1:360
    q=quiver3(x(inds==i),y(inds==i),z(inds==i),UX1(inds==i),UY1(inds==i),UZ1(inds==i),.07,'linestyle','--','color',cmap(i,:),'linewidth',1.5); %32k
    q.ShowArrowHead = 'off';
end
view(-90,0)



