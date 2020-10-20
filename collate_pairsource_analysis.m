%% Collate pair source analysis data
%
%   uses pair source analysis data from collision_atom_distribution.m
%
% DKS 2020

% NOTE: mJ indices run 0,1 for pair source (but 1,0,-1 for later data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% variables to collate
% * total atom number, total scattered number
%
% * number statistics: hemisphere, single-mode localised 
%
% * g2
%
% * mode analysis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%% Configure
% path_data   = 'C:\Users\David\Dropbox\PhD\diagrams\data_analysis\k_space\full_analysis\20200612';
path_data   = '/home/david/Documents/collision/pairsrc_out';


file_head   = 'anal_run';


%% collate data
S_all = load_collate(path_data,file_head,1);

S = cell2mat(S_all);
n_data = length(S);

clearvars S_all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% collate vars and process
%% Total atom number, total scattered number
Nsc_mean = [S.N_sc_mean];
Nsc_std = [S.N_sc_std];
Ntot_mean = [S.N_tot_mean];
Ntot_std = [S.N_tot_std];

%%% fit -- only a guide
vm=[Ntot_mean.^2]';      %vandermode matrix
beta = vm\Nsc_mean';

xx_NN = linspace(0,3.5e3);
yy_NN = beta*xx_NN.^2;


%% correlation g2
corr = [S.corr];

g2_amp = [corr.g2_amp];
g2_amp_err = [corr.g2_amp_err];

sig_bb = [corr.sig_bb];
sig_bb_err = [corr.sig_bb_err];


%% Halo characteristics
halo = [S.halo];

w = [halo.width];
w_err = [halo.width_err];

n_mocc = [halo.n_mocc];
n_mocc_err = [halo.n_mocc_err];

M = [halo.M];
M_err = [halo.M_err];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA VIS
%% Tuning total scattered number via condensate N
h_NN = figure('name','tune_scattered_number');
p_NN = ploterr(Ntot_mean,Nsc_mean,Ntot_std,Nsc_std,'ro','hhxy',0);
p_NN(1).MarkerFaceColor='w';

hold on;
p_NN_line = plot(xx_NN,yy_NN,'k-');
uistack(p_NN_line,'bottom');

xlabel('total detected counts');
ylabel('scattered counts');
box on;


%% PLOT -- halo mode characteristic vs. scat num
h_mode = figure('Name','mode');

% halo rms width
subplot(3,1,1);
hold on;
p_width = ploterr(Nsc_mean,w,Nsc_std,w_err,'ks','hhxy',0);
p_width(1).MarkerFaceColor = p_width(2).Color;
ylabel('$\delta k$');
box on;

% BB correlation length
subplot(3,1,2);
hold on;
p_sigbb = ploterr(Nsc_mean,sig_bb,Nsc_std,sig_bb_err,'ko','hhxy',0);
p_sigbb(1).MarkerFaceColor = 'w';
ylabel('$\sigma_\textrm{BB}$');
box on;

% Number of modes
subplot(3,1,3);
hold on;
p_M = ploterr(Nsc_mean,M,Nsc_std,M_err,'k^','hhxy',0);
p_M(1).MarkerFaceColor = 'w';

xlabel('scattered counts');
ylabel('$M$');
box on;

%% PLOT -- g2 vs mode occupancy
h_g2n = figure('Name','g2_mocc');
hold on;

col_data = 'k';
col_biv = 0.8*ones(1,3);


% data
p_g2n = ploterr(n_mocc,g2_amp-1,n_mocc_err,g2_amp_err,'o','hhxy',0);
p_g2n(1).MarkerFaceColor=col_data;
set(p_g2n,'Color',col_data);

ax=gca;
ax.XScale = 'log';
ax.YScale = 'log';

% theory
g2n.xx = geospace_lim(ax.XLim,1e2);
g2n.yy = 2 + 1./g2n.xx;
g2n.p = plot(g2n.xx,g2n.yy - 2,'-','Color',col_data,'DisplayName','theory');
uistack(g2n.p,'bottom');

xlabel('mode occupancy');
ylabel('correlation amplitude, $g^{(2)}_{\uparrow\downarrow} - 2$');
box on;

% BIV shaded region
g2_biv = 3 + 2*sqrt(2);

y0 = ax.YLim(1);
y1 = g2_biv-2;
p_biv = patch([ax.XLim fliplr(ax.XLim)],[y0,y0,y1,y1],col_biv,'EdgeColor','none');
uistack(p_biv,'bottom');

ax.Layer = 'top';

%% save
out_str = ['collate_',file_head,'_',getdatetimestr];
save(fullfile(path_data,out_str));