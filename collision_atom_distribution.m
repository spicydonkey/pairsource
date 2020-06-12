%% Characterising k-space atom number distribution
% DKS 2020
%
%
% compare:
%   * mJ
%   * source and theta/2 rotation


%% config
%%% DATA
% data_path = 'C:\Users\David\data\2020_revisit\kspace_hotspot_filt\out\collate_collate_config_20200603_165608.mat';
% data_path = 'C:\Users\David\data\bell_momentumspin\pairsource\90\rev\run_1_rev.mat';
% b_flip_mJ = true;

% data structure contextual
idx_parT_anal = 1;        % FOR BELL YROT: 1 = pi/2; 2 = source

%%% 1. Radial dist
n_bin_r = 200;
r_lim = [0.8,1.2];

%%% 2. Az-el spherical distribution
n_bin_azel =[100,50];

%%% 3. Number statistics (hemisphere)
hemi.theta = pi/2;
hemi.az = pi;
hemi.el = 0;

%%% 4. k-localised number statistics
kzone.theta = 0.05;     % 0.03*2^0;

%%% 5. Quantum efficiency
zone.theta = pi/15;     % cone half-angle

%%% 6. g2 correlation
corr.nbin_th = 5e2;


%%% STAT
bs.nboot = 100;
bs.nboot_g2 = 25;
bs.frac_g2 = 0.2;


%%% params
z_trunc = 0.8;      % may be overridden by config
det_qe = 0.1;


%%% Graphics and annotation
col_mJ = {'b','r'};
str_mJ = {'$1\quad|\!\uparrow\rangle$','$0\quad|\!\downarrow\rangle$'};
mark_mJ = {'^','v'};
str_g2 = {'$\uparrow\uparrow$','$\uparrow\downarrow$';'$\downarrow\uparrow$','$\downarrow\downarrow$'};
col_g2 = {'b','k';'k','r'};

% check mJ indexing
if ~exist('b_flip_mJ','var')
    b_flip_mJ = false;
end
if b_flip_mJ
    col_mJ = fliplr(col_mJ);
    str_mJ = fliplr(str_mJ);
    mark_mJ = fliplr(mark_mJ);
    str_g2 = fliplr(flipud(str_g2));
    col_g2 = fliplr(flipud(col_g2));
    idx_mJ0 = 1;
else
    idx_mJ0 = 2;
end


%% Preprocess
S = load(data_path);                % load data

K = S.k_par{idx_parT_anal};     

% %%% DEBUG START
% K = K(1:round(size(K,1)/3),:);
% disp('DEBUGGING: reducing data for speed.');
% %%% DEBUG END

n_shot = size(K,1);
n_mJ = size(K,2);                   % # of internal states
nk = shotSize(K);
nk_max = max(nk(:));                % max count in a mJ halo

K_coll = cell_vertcat(K);       


%% 1. radial distribution
%%% Configure radial histogram.
r_ed=linspace_lim(r_lim,n_bin_r+1);
r_cent = edge2cent(r_ed);

%%% Bootstrapping way.
k_norm = cellfun(@(x) vecnorm(x,2,2),K,'uni',0);

fun_nr = @(x) rdist_truncsph(vertcat(x{:}),r_ed,z_trunc)/numel(x);

n_r = arrayfun(@(I) fun_nr(k_norm(:,I)),1:2,'uni',0);
[n_r_err,~,n_r_bs] = arrayfun(@(I) bootstrap(bs.nboot,fun_nr,k_norm(:,I)),1:2,'uni',0);

%%% FIT radial dist
kdist_fit = cell(1,2);
kdist_mdl = cell(1,2);

k_width = NaN(1,2);
k_width_fiterr = NaN(1,2);
k_width_bs = cell(1,2);
k_width_staterr = NaN(1,2);

p_const = cell(1,4);        % all params unconstrained
for ii = 1:2
    ty = n_r{ii};
    p0 = [max(ty),1,0.03,1];
    [~,kdist_fit{ii},kdist_mdl{ii}] = fit_gauss_1d(r_cent,ty,p0,p_const,0);
    
    % get halo width
    [k_width(ii),k_width_fiterr(ii)] = get_gauss1d_width(kdist_mdl{ii});
    % bootstrapped error
    k_width_bs{ii} = arrayfun(@(I) get_gauss1d_width(mdl_fit_gauss_1d(...
        r_cent,n_r_bs{ii}(I,:),p0,p_const,0)),1:bs.nboot);
    k_width_staterr(ii) = std(k_width_bs{ii});
end
k_width_err = quadsum(k_width_staterr,k_width_fiterr);


%% PLOT radial distribution
h_rdist = figure('name','rdist');
hold on;

% p_rdist = gobjects(1,2);
p_rdist = gobjects(2,2);
p_rdist_fit = gobjects(1,2);
for ii = 1:2
    p_rdist(ii,:) = ploterr(r_cent,n_r{ii},[],n_r_err{ii},col_mJ{ii});
    set(p_rdist(ii,1),{'LineStyle','Marker','MarkerFaceColor','DisplayName'},{'none','.','w',str_mJ{ii}});
    
    % fit
    p_rdist_fit(ii) = plot(kdist_fit{ii}.x,kdist_fit{ii}.y,col_mJ{ii});
    uistack(p_rdist_fit(ii),'bottom');
end
xlabel('momentum k');
ylabel('count density (arb. u.)');

ax=gca;
ax.Layer = 'top';
grid on;

lgd = legend([p_rdist([1,2])]);
lgd.Title.String = '$m_J$';
box on;


%% 2. spherical density distribution
%%% setup analysis params
az_ed = linspace(-pi,pi,n_bin_azel(1)+1);
el_ed = linspace(-pi/2,pi/2,n_bin_azel(2)+1);
az_c = edge2cent(az_ed);
el_c = edge2cent(el_ed);

% grid
[AZ_ed,EL_ed] = ndgrid(az_ed,el_ed);
[AZ_c,EL_c] = ndgrid(az_c,el_c);
b_trunc = abs(EL_c) > asin(z_trunc);


%%% histogram
norm_azel = 1./(n_shot*solang_azel_grid(AZ_ed,EL_ed));      % normalisation
n_azel_orig = cellfun(@(x) norm_azel.*hist_zxy_azel(x,{az_ed,el_ed}),K_coll,'uni',0);

% truncate region to NaN
n_azel = n_azel_orig;
for ii = 1:numel(n_azel)
    n_azel{ii}(b_trunc) = nan;
end
n_azel_se = cellfun(@(x) sqrt(x)/sqrt(n_shot),n_azel,'uni',0);


% averages in azim/polar axes
% TODO -- what is the statistical uncertainty from model
n_az_avg = cellfun(@(x) mean(x,2,'omitnan'),n_azel,'uni',0);
n_az_std = cellfun(@(x) std(x,[],2,'omitnan'),n_azel,'uni',0);

n_el_avg = cellfun(@(x) mean(x,1,'omitnan'),n_azel,'uni',0);
n_el_std = cellfun(@(x) std(x,[],1,'omitnan'),n_azel,'uni',0);



%%% PLOT - 2D AZEL dist
h_azel = gobjects(1,2);
p_azel = gobjects(1,2);
cb = gobjects(1,2);
for ii = 1:2
    h_azel(ii) = figure('name',['azel_',num2str(ii)]);
    p_azel(ii) = imagesc(az_c,el_c,n_azel{ii}');
    ax = gca;
    ax.YDir = 'normal';
    
    xlabel('$\theta$');
    ylabel('$\phi$');
    
    cb(ii) = mycolorbar;
    cb(ii).Label.String = 'count density (sr$^{-1}$)';
    
    axis equal tight;
end


%% 2.1. Azimuthal distribution
%%% PLOT
h_azim = figure('name','azim');
hold on;
p_azim = cell(1,2);
for ii = 1:2 
    p_azim{ii} = shadedErrorBar(az_c,n_az_avg{ii},n_az_std{ii}/sqrt(n_bin_azel(2)),...
        'lineprops',col_mJ{ii});
end
p_azim = cell2mat(p_azim);

xlabel('azimuthal angle $\theta$');
ylabel('count density (sr$^{-1}$)');

ax=gca;
xlim(pi*[-1,1]);
ax.YLim(1) = 0;


%% 2.2. Polar distribution
%%% PLOT
h_elev = figure('name','elev');
hold on;
p_elev = cell(1,2);
for ii = 1:2 
%     p_elev(ii) = plot(vel,nk_azel{ii}(idx_az_plot,:));
%     p_elev(ii) = plot(el_c,n_azel{ii}(idx_az_plot,:));
    p_elev{ii} = shadedErrorBar(el_c,n_el_avg{ii},n_el_std{ii}/sqrt(n_bin_azel(1)),'lineprops',col_mJ{ii});
end
p_elev = cell2mat(p_elev);

xlabel('elevation angle $\phi$');
ylabel('count density (sr$^{-1}$)');

ax = gca;
ax.YLim(1) = 0;
xlim(0.5*pi*[-1,1]);


%% 3. number statistics -- halo-hemisphere integrated
% NOTE: counts must be in hemisphere to sample from just one part of correlated pair
% get counts in azel grid
[~,~,nn] = cellfun(@(x) inCone(x,hemi.az,hemi.el,hemi.theta),K,'uni',0);
hemi.nk = cell2mat(nn);

hemi.n_max = maxall(hemi.nk);           % max counts in hemisphere
hemi.n_ed = -0.5:(hemi.n_max+1);
hemi.n_ct = edge2cent(hemi.n_ed);

% statistics from data
hemi.nk_avg = mean(hemi.nk);
hemi.nk_std = std(hemi.nk);

% histogram
fun_countpdf = @(x) histcounts(x,hemi.n_ed,'Normalization','pdf');
hemi.p_nk = arrayfun(@(I) fun_countpdf(hemi.nk(:,I)),[1,2],'uni',0);
hemi.p_nk_err = arrayfun(@(I) bootstrap(bs.nboot,fun_countpdf,hemi.nk(:,I)),...
    [1,2],'uni',0);


%% PLOT time series -- experiment stability
h_ntseries = figure('name','n_tseries');
hold on;
p_hemi_nt = plot(hemi.nk);
xlabel('shot number');
ylabel('atom counts in hemisphere');


%% PLOT number statistics
h_nstat = figure('name','n_stat');
hold on;

p_nstat = gobjects(2,2);
for ii = 1:2
    p_nstat(ii,:) = ploterr(hemi.n_ct,hemi.p_nk{ii},[],hemi.p_nk_err{ii},col_mJ{ii});
    set(p_nstat(ii,1),{'Marker','LineStyle','DisplayName'},...        
        {'.','none',str_mJ{ii}});
end

% annotate
ax=gca;
ax.Layer = 'top';

xlabel('atom counts in halo');
ylabel('probability');

lgd = legend(p_nstat([1,2]),'location','best');


%% Time segment
% % hemi.idx_shot_seg = 1.3e4 + (0:2000);
% % hemi.idx_shot_seg = [1:700,750:n_shot];
% hemi.idx_shot_seg = 1:n_shot;        % DO ALL
% if length(hemi.idx_shot_seg)<n_shot
%     disp('Analysing a segment of data.');
% end
%     
% hemi.nk = hemi.nk(hemi.idx_shot_seg,:);
% 
% % histogram
% hemi.N_nk = arrayfun(@(I) nhist(hemi.nk(:,I),{hemi.n_ed}),1:2,'uni',0);
% hemi.p_nk = cellfun(@(x) x./sum(x),hemi.N_nk,'uni',0);    % normalise
% 
% 
% %%% PLOT timeseries
% h_ntseries_seg = figure('name','n_tseries_seg');
% plot(hemi.idx_shot_seg,hemi.nk');
% xlabel('shot number');
% ylabel('atom counts in halo');
% 
% 
% %%% PLOT counting statistics
% figure(h_nstat);
% hold on;
% 
% p_nstat_seg = gobjects(1,2);
% for ii = 1:2
% %     p_nstat_seg(ii) = stairs(n_ct,p_nk{ii},'Color',col{ii},...
% %         'LineStyle','--','DisplayName','segment');
%     p_nstat_seg(ii) = plot(hemi.n_ct,hemi.p_nk{ii},'Marker',mark_mJ{ii},'Color',col_mJ{ii},...
%         'LineStyle','none','DisplayName','segment');
% end


%% Theoretical distributions
%%% Thermal
thermal.xx = hemi.n_ct;
thermal.yy = arrayfun(@(mu) thermalpdf(hemi.n_ct,mu),hemi.nk_avg,'uni',0);

%%% Poissonian 
poiss.xx = hemi.n_ct;
poiss.yy = arrayfun(@(L) poisspdf(hemi.n_ct,L),hemi.nk_avg,'uni',0);

%%% Multimode thermal 
fun_mmthermal{1} = @(M,x) nbinpdf_mmthermal(x,hemi.nk_avg(1),M);
fun_mmthermal{2} = @(M,x) nbinpdf_mmthermal(x,hemi.nk_avg(2),M);

b_mmthermal_ok = true;
mmthermal.fitmdl=cell(1,2);
for ii = 1:2
    try
        mmthermal.fitmdl{ii} = fitnlm(hemi.n_ct,hemi.p_nk(:,ii),fun_mmthermal{ii},10);
    catch
        b_mmthermal_ok = false;
        disp('mmthermal fit has failed.');
    end
end
mmthermal.xx = linspace(0,max(hemi.n_ct),1e3);
if b_mmthermal_ok
    mmthermal.yy = cellfun(@(x) feval(x,mmthermal.xx),mmthermal.fitmdl,'uni',0);
else
    mmthermal.yy = cellfun(@(x) NaN(size(mmthermal.xx)),mmthermal.fitmdl,'uni',0);
end


%% PLOT theoretical distributions
figure(h_nstat);
hold on;

p_nstat_therm = gobjects(1,2);
p_nstat_poiss = gobjects(1,2);
p_nstat_mmthermal = gobjects(1,2);

for ii = 1:2
    p_nstat_therm(ii) = plot(thermal.xx,thermal.yy{ii},...
        'Color',col_mJ{ii},'LineStyle',':','DisplayName','thermal');
    p_nstat_poiss(ii) = plot(poiss.xx,poiss.yy{ii},...
        'Color',col_mJ{ii},'LineStyle','--','DisplayName','Poisson');
    p_nstat_mmthermal(ii) = plot(mmthermal.xx,mmthermal.yy{ii},'Color',col_mJ{ii},...
        'linestyle','-','DisplayName','MM thermal');
    
    % layer order
    uistack(p_nstat_therm(ii),'bottom');
    uistack(p_nstat_poiss(ii),'bottom');
    uistack(p_nstat_mmthermal(ii),'bottom');
end

set(p_nstat,'MarkerFaceColor','w');
set(gca,'YScale','log');
box on;


%% 4. Number statistics -- k-localised
%%% configure
% zone params
kzone.vaz = 0:2*kzone.theta:2*pi;
kzone.vel = 0;


%%% setup
[kzone.az,kzone.el] = ndgrid(kzone.vaz,kzone.vel);
kzone.n_zone = numel(kzone.az);


%%% check if any overlap in zones
dth_pair_matrix = NaN([size(kzone.az),size(kzone.az)]);
for ii = 1:numel(dth_pair_matrix)
    [ii1,ii2,ii3,ii4] = ind2sub(size(dth_pair_matrix),ii);
    if (ii1==ii3) && (ii2==ii4)
        continue
    end
    dth_pair_matrix(ii1,ii2,ii3,ii4) = diffAngleSph(kzone.vaz(ii1),kzone.vel(ii2),kzone.vaz(ii3),kzone.vel(ii4));
end
nzone_overlap = sumall(dth_pair_matrix<(1.99*kzone.theta))/2;
if nzone_overlap>0
    warning(sprintf('there are %d overlapping zones',nzone_overlap));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% get number in each k-zone
kzone.n_k = cell(size(kzone.az));
for ii = 1:kzone.n_zone
    taz = kzone.az(ii);
    tel = kzone.el(ii);
    [~,~,nn] = cellfun(@(x) inCone(x,taz,tel,kzone.theta),K,'uni',0);
    kzone.n_k{ii} = cell2mat(nn);
end

% statistics
kzone.n_k_avg = cell2mat(cellfun(@(x) mean(x),kzone.n_k,'uni',0));
kzone.n_k_std = cell2mat(cellfun(@(x) std(x),kzone.n_k,'uni',0));

kzone.n_ed = -0.5:(maxall(cell2mat(kzone.n_k))+1);
kzone.n_ct = edge2cent(kzone.n_ed);

%%% histogram
kzone.n_nk = cell(1,2);
% kzone.p_nk = cell(1,2);
for ii = 1:2
    kzone.n_nk{ii} = nan([size(kzone.az),length(kzone.n_ct)]);
    for jj = 1:length(kzone.vaz)
        for kk = 1:length(kzone.vel)
            kzone.n_nk{ii}(jj,kk,:) = nhist(kzone.n_k{jj,kk}(:,ii),{kzone.n_ed});
        end
    end
end
%normalise to give probability
kzone.p_nk = cellfun(@(x) x./sum(x,3),kzone.n_nk,'uni',0);

% statistics
kzone.p_nk_avg = cellfun(@(x) squeeze(mean(x,[1,2])),kzone.p_nk,'uni',0);
kzone.p_nk_se = cellfun(@(x) squeeze(std(x,[],[1,2]))./sqrt(kzone.n_zone),kzone.p_nk,'uni',0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Plot k-localised counting stats
h_nstat_kloc = figure('name','nstat_kloc');
hold on;

p_kloc = cell(1,2);
for ii = 1:2
    p_kloc{ii} = ploterr(kzone.n_ct',kzone.p_nk_avg{ii},[],kzone.p_nk_se{ii},...
        [col_mJ{ii},'o']);
    set(p_kloc{ii}(1),{'MarkerFaceColor','DisplayName'},{'w',str_mJ{ii}});
end

ax=gca;
ax.Layer = 'top';
ax.YScale = 'log';
xlabel('count in zone');
ylabel('probability');
grid on;

lgd = legend([p_kloc{1}(1),p_kloc{2}(1)]);
lgd.Title.String = '$m_J$';
title(sprintf('$\\alpha =$ %0.3f,\\quad%d zones',kzone.theta,kzone.n_zone));


%% 5. Quantum efficiency
%%% configure
zone.vaz = 0:2*zone.theta:(2*pi-0.01);    % there will be overlap!
zone.vaz = zone.vaz(1:end-1);
zone.vel = 0;

% get opposite zone
zone.vaz_bb = wrapToPi(zone.vaz+pi);
zone.vel_bb = -zone.vel;

% define reference zone
[AZ_zone,EL_zone] = ndgrid(zone.vaz,zone.vel);
n_azel_zone = numel(AZ_zone);

% opposite zones
[AZ_zone_bb, EL_zone_bb] = ndgrid(zone.vaz_bb,zone.vel_bb);
 
zone.az = {AZ_zone, AZ_zone_bb};
zone.el = {EL_zone, EL_zone_bb};

% Get overlapping zones.
zone.dth = diffAngleSph(zone.az{1},zone.el{1},zone.az{2}',zone.el{2}');
zone.b_overlap = abs(zone.dth) < 2*zone.theta;      % there is overlap
zone.b_diag = logical(eye(n_azel_zone));            % will be exactly BB

%%% Count atoms in zones
zone.Nm = cell(1,2);
for ii = 1:2
    zone.Nm{ii} = cell(size(AZ_zone));       
    for jj = 1:n_azel_zone
        taz = zone.az{ii}(jj);
        tel = zone.el{ii}(jj);
        [~,~,nn] = cellfun(@(x) inCone(x,taz,tel,zone.theta),K,'uni',0);
        nn = cell2mat(nn);
        zone.Nm{ii}{jj} = nn;
    end
end

% sum over mJ in k-zone
zone.N = arrayfun(@(I) cellfun(@(x) sum(x,2),zone.Nm{I},'uni',0),1:2,'uni',0);
% zone.Navg = arrayfun(@(I) cellfun(@(x)  mean(x),zone.N{I}),1:2,'uni',0);

% Evaluate normalised variance of number difference.
zone.varnorm = NaN(n_azel_zone);
zone.varnorm_err = NaN(n_azel_zone);
zone.varnorm_bs_mean = NaN(n_azel_zone);
for ii = 1:n_azel_zone
    for jj = 1:n_azel_zone
        tN = packrows([zone.N{1}{ii}, zone.N{2}{jj}]);
        
        % All data
        zone.varnorm(ii,jj) = normDiffVar(tN);
        % Bootstrap to find error of NDV.
        [zone.varnorm_err(ii,jj),zone.varnorm_bs_mean(ii,jj)] = bootstrap(bs.nboot,@normDiffVar,tN);
    end
end

% % evaluate number difference between pairs of zones
% zone.Ndiff = cell(n_azel_zone,n_azel_zone);
% for ii = 1:n_azel_zone
%     for jj = 1:n_azel_zone
%         zone.Ndiff{ii,jj} = zone.N{1}{ii} - zone.N{2}{jj};
%     end
% end
% 
% %%% evaluate statistics
% zone.Ndiff_avg = cellfun(@(x) mean(x),zone.Ndiff);
% zone.Navg_sum = zone.Navg{1} + zone.Navg{2}';
% 
% zone.Ndiff_var = cellfun(@(x) var(x),zone.Ndiff);
% % zone.Ndiffsq_var = cellfun(@(x) var(x.^2),zone.Ndiff);
% 
% zone.varnorm = zone.Ndiff_var./zone.Navg_sum;
% % zone.varnorm_err = sqrt(zone.Ndiffsq_var./(zone.Navg_sum*n_shot));        % error est from poisson assumption

% set same-bin analysis to nan
zone.varnorm(zone.b_overlap) = NaN;
zone.varnorm_err(zone.b_overlap) = NaN;

%%% stat
zone.varnorm_corr_mean = mean(zone.varnorm(zone.b_diag),'omitnan');
zone.varnorm_unc_mean = mean(zone.varnorm(~zone.b_diag),'omitnan');
zone.varnorm_corr_serr = std(zone.varnorm(zone.b_diag),'omitnan')/sqrt(sumall(zone.b_diag));
zone.varnorm_unc_serr = std(zone.varnorm(~zone.b_diag),'omitnan')/sqrt(sumall(~zone.b_diag));


%%% Estimate detector QE
qe_meas = 1 - zone.varnorm_corr_mean;
qe_meas_err = zone.varnorm_corr_serr;
qe_meas_ferr = qe_meas_err/qe_meas;


%% PLOT norm variance
h_varnorm = figure('name','varnorm');
hold on;
p_corr = ploterr(zone.vaz,zone.varnorm(zone.b_diag),[],zone.varnorm_err(eye(n_azel_zone)==1),'ko','hhxy',0);
set(p_corr(1),'MarkerFaceColor','k');
p_unc = ploterr(0,zone.varnorm(~zone.b_diag),[],zone.varnorm_err(eye(n_azel_zone)~=1),'ko','hhxy',0);
set(p_unc(1),'MarkerFaceColor','w');

% mean
ax=gca;
shadedErrorBar(ax.XLim,zone.varnorm_corr_mean*[1,1],zone.varnorm_corr_serr*[1,1],...
    'lineprops','g');
shadedErrorBar(ax.XLim,zone.varnorm_unc_mean*[1,1],zone.varnorm_unc_serr*[1,1],...
    'lineprops','k');

%%% annotate
set(p_unc,'Visible','off');
uistack(p_unc,'top');
uistack(p_corr,'top');

ax.Layer = 'top';
xlabel('zone azimuthal angle $\theta$');
ylabel('normalised variance');


%% 2D plot
h_varnorm2d = figure('name','varnorm_2d');

imagesc(zone.varnorm,'AlphaData',double(~isnan(zone.varnorm)));

ax=gca;
axis equal tight;
ax.XAxisLocation='top';
xlabel('BB zone $j$');
ylabel('zone $i$');
cb = mycolorbar;
cb.Label.String = 'normalised variance';
colormap(flipud(viridis));


%% 3D bar
% h_varnorm3d = figure('name','varnorm_3d');
% b_vecnorm = bar3(zone.varnorm);
% 
% for ii = 1:length(b_vecnorm)
%     b_vecnorm(ii).CData = b_vecnorm(ii).ZData;
% end
% set(b_vecnorm,'FaceColor','interp');
% cb=mycolorbar;



%% 6. g2 analysis
%%% configure pairwise difference angle histogram
corr.th_ed = linspace(0,pi,corr.nbin_th+1);
corr.th_ct = edge2cent(corr.th_ed);


%%% angular g2
fun_g2 = @(x) angular_g2(x,corr.th_ed);     % Functionalised for bootstrap

corr.g2 = cell(2);
corr.g2_err = cell(2);
corr.g2_bs_mean = cell(2);
corr.g2_bs = cell(2);           % find fit amp error by BS
for ii = 1:2
    for jj = 1:2
        tK = packrows(K(:,unique([ii,jj])));            % pack to group mJ types
        corr.g2{ii,jj} = fun_g2(tK);                    % super simple!
        [corr.g2_err{ii,jj},corr.g2_bs_mean{ii,jj},corr.g2_bs{ii,jj}] = bootstrap(bs.nboot_g2,fun_g2,tK,bs.frac_g2);
    end
end

%% Fit model
% models
fun_g2_bb = @(b,x) 1 + b(1)*exp(-(x-pi).^2/(2*b(2).^2));
fun_g2_cl = @(b,x) 1 + b(1)*exp(-x.^2/(2*b(2).^2));
fun_g2_bb_cl = @(b,x) 1 + b(1)*exp(-(x-pi).^2/(2*b(2).^2)) + b(3)*exp(-x.^2/(2*b(4).^2));

% initial params
sigk0 = 0.03;           
par0_cl = [1,sigk0];

% fit
corr.g2_mdl = cell(2);
for ii = 1:2
    for jj = 1:2
        if ii == jj
            corr.g2_mdl{ii,jj} = fitnlm(corr.th_ct,corr.g2{ii,jj},fun_g2_cl,par0_cl);
        else
            tp0 = [corr.g2{ii,jj}(end),sigk0];
            corr.g2_mdl{ii,jj} = fitnlm(corr.th_ct,corr.g2{ii,jj},fun_g2_bb,tp0);
        end
    end
end
% evaluate fit
corr.th_fit = linspace(0,pi,1e4);
corr.g2_fit = cellfun(@(m) feval(m,corr.th_fit),corr.g2_mdl,'uni',0);

%%% get g2_BB params
% NOTE: fit error from all is ~0.1 of BS error
% From all data.
[corr.g2_amp, corr.sig_bb, corr.g2_amp_fiterr, corr.sig_bb_fiterr] = get_g2bb_fitparams(corr.g2_mdl{1,2});

% Bootstrapping fit errors.
tg2_bs = corr.g2_bs{1,2};
tp0_bs = [corr.g2_amp,corr.sig_bb];
[corr.g2_amp_bs, corr.sig_bb_bs] = arrayfun(@(I) get_g2bb_fitparams(...
    fitnlm(corr.th_ct,tg2_bs(I,:),fun_g2_bb,tp0_bs)),1:bs.nboot_g2);
corr.g2_amp_staterr = sqrt(bs.frac_g2)*std(corr.g2_amp_bs);
corr.sig_bb_staterr = sqrt(bs.frac_g2)*std(corr.sig_bb_bs);

% effective error
corr.g2_amp_err = quadsum(corr.g2_amp_staterr,corr.g2_amp_fiterr);
corr.sig_bb_err = quadsum(corr.sig_bb_staterr,corr.sig_bb_fiterr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT
h_g2 = figure('name','g2');
hold on;

p_g2 = gobjects(2,2,2);
p_fit = gobjects(2,2);
for ii = 1:2
    for jj = 1:2
        % data and fit
        p_g2(ii,jj,:) = ploterr(corr.th_ct,corr.g2{ii,jj},[],corr.g2_err{ii,jj},[col_g2{ii,jj},'.']);
        p_fit(ii,jj) = plot(corr.th_fit,corr.g2_fit{ii,jj},col_g2{ii,jj});
        
        % annotate
        set(p_g2(ii,jj,:),{'DisplayName'},{str_g2{ii,jj}});
    end
end
ax = gca;
ax.Layer='top';
xlabel('relative angle $\theta$');
ylabel('$g^{(2)}$');
xlim([0,pi]);
box on;
[~,data_name]=fileparts(data_path);
title(data_name,'Interpreter','none');

set(p_g2(2,1,:),'Visible','off');
set(p_fit(2,1),'Visible','off');

lgd = legend([p_g2(1,1),p_g2(2,2),p_g2(1,2)],'Location','best');


%% 7. BEC vs scattered number
if isfield(S,'txy')
    zxy = txy2zxy(S.txy);
    
    %%% total number
    % simplest approximation for number in single condensate
    N_tot = shotSize(zxy);
    N_tot_mean = mean(N_tot);        % heavily saturated
    N_tot_std = std(N_tot);
end
%%% total scattered number
ratio_trunc = z_trunc;      % truncation ratio is just z_trunc

% counts in truncation volume
N_sc = sum(nk,2);
N_sc_mean = mean(N_sc);
N_sc_std = std(N_sc);
N_sc_ferr = (N_sc_std/sqrt(n_shot))/N_sc_mean;      % stderr in mean


%% 8. Mode analysis
% REF:  SOM, Hodgman, et al. PRL 2017.
halo.sig_k = corr.sig_bb/1.1;
halo.sig_k_err = corr.sig_bb_err/1.1;      
halo.width = k_width(idx_mJ0);
halo.width_err = k_width_err(idx_mJ0);          

halo.Vs = 4*pi*sqrt(2*pi)*halo.width;       % approx vol of 3D spherical shell with gaussian r-dist
halo.Vm = (2*pi)^(3/2)*halo.sig_k^3;        % coherence vol (3D gaussian)

halo.M = halo.Vs/halo.Vm;      % modes PER halo (x2 for both halos)
halo.M_ferr = halo.width_err/halo.width + 3*halo.sig_k_err/halo.sig_k;   % fractional errors
halo.M_err = halo.M * halo.M_ferr;

% mode occupation
halo.n_mocc = N_sc_mean/(2*halo.M*ratio_trunc*det_qe);
halo.n_mocc_ferr = N_sc_ferr + halo.M_ferr;
halo.n_mocc_err = halo.n_mocc * halo.n_mocc_ferr;


%% save data
% cull large chunk of rawish data
clearvars S zxy;

% save
out_str = ['anal_',data_name,'_',getdatetimestr];
save(out_str);


%% Local functions
%%% radial dist
function mdl = mdl_fit_gauss_1d(x,y,p0,pConst,verbose)
    [~,~,mdl] = fit_gauss_1d(x,y,p0,pConst,verbose);
end

function [w,w_err] = get_gauss1d_width(mdl)
    w = mdl.Coefficients.Estimate(3);
    w_err = mdl.Coefficients.SE(3);
end


%%% QE
% normalised difference variance.
%   x is a Nx1 cell array of 1x2 double vector of correlated signal.
function v = normDiffVar(x)
x = unpackrows(x);
dx = x(:,2) - x(:,1);       
x_avg = mean(x,1);
v = var(dx)/sum(x_avg);
end


%%% g2
function g2 = angular_g2(K,th_ed)
K = unpackrows(K);
[N_dth_corr,N_dth_all] = pairAngCount(K,th_ed,1);
% normalise
nn_dth_corr = N_dth_corr./sumall(N_dth_corr);
nn_dth_all = N_dth_all./sumall(N_dth_all);
g2 = nn_dth_corr./nn_dth_all;
end

function [amp,sig,amp_se,sig_se] = get_g2bb_fitparams(mdl)
    amp = mdl.Coefficients.Estimate(1);
    sig = mdl.Coefficients.Estimate(2);
    amp_se = mdl.Coefficients.SE(1);
    sig_se = mdl.Coefficients.SE(2);
end
