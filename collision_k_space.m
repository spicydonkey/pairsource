%% Collision schematic in k-space
% 
% DKS 2020

%% config
shot_single = 14000;
shot_collate = shot_single+(1:100);      


%% preproc data
% data_preproc = 'C:\Users\David\Dropbox\PhD\diagrams\data_analysis\k_space\exp1_collated_20180815_1.mat';
data_preproc = 'C:\Users\David\data\2020_revisit\kspace_hotspot_filt\config_1_collate_20200602_170905.mat';
S = load(data_preproc);


%% get data
K = S.k_par{1};       % get k-space data


%% plot
h_k3d = figure('name',['k3d_',num2str(shot_single)],'Renderer','painters');

p_collate = plot_zxy(K(shot_collate,:));
p_single = plot_zxy(K(shot_single,:));

%%% graphics
mark_size = 30;
mark_type = 'o';
coll_alpha = 0.2;               % transparency for collate data
single_alpha = 0.8;
col0_rgb = [0,0,1; 1,0,0];      % orig blue/red
col0_hsl = rgb2hsl(col0_rgb);     
coll_hsl = col0_hsl;            % lighter shade
coll_hsl(:,3) = 0.8;
coll_rgb = hsl2rgb(coll_hsl);

% collate data to unsaturated and transparent
set(p_collate,'marker',mark_type);
set(p_collate,'sizedata',mark_size);
set(p_collate,'MarkerEdgeColor','none');
set(p_collate(1),'MarkerFaceColor',coll_rgb(1,:));
set(p_collate(2),'MarkerFaceColor',coll_rgb(2,:));
set(p_collate,'MarkerFaceAlpha',coll_alpha);

% single shot
set(p_single,'marker',mark_type);
set(p_single,'sizedata',mark_size);
set(p_single,'MarkerEdgeColor','none');
set(p_single(1),'MarkerFaceColor',col0_rgb(1,:));
set(p_single(2),'MarkerFaceColor',col0_rgb(2,:));
set(p_single,'MarkerFaceAlpha',single_alpha);

%%% sphere
% SURF way
% n_s = 30;       % num of sphere face
% [Xs,Ys,Zs] = sphere(n_s);
% Ss = surf(Xs,Ys,Zs,'FaceColor','none','EdgeAlpha',0.5);

% SCAFFOLD way
ns_sparse = 6;
ns_dense = 3*ns_sparse;
l_sph_sparse = plot_scaffold(sphere_scaffold(ns_sparse));
set(l_sph_sparse,'Color',[0 0 0 0.5]);
l_sph_dense = plot_scaffold(sphere_scaffold(ns_dense));
set(l_sph_dense,'Color',[0 0 0 0.1]);


%%% cartesian axes from origin
ax_lim = 1.5;
lcart = gobjects(1,3);
lcart(1) = line([0,ax_lim],[0,0],[0,0]);
lcart(2) = line([0,0],[0,ax_lim],[0,0]);
lcart(3) = line([0,0],[0,0],[0,ax_lim]);
set(lcart,'Color','k');


%%% annotate
ax=gca;
xlabel('$k_x$');
ylabel('$k_y$');
zlabel('$k_z$');
set(ax, 'Projection','perspective');
view([142.5,24]);
