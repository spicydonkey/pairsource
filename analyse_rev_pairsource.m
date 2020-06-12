% reanalyse g2 in angular coord
% DKS 2020

clear all;

tic;

%% config
% data_fname = {'run_10_rev.mat'};
% data_fname = {'run_1_rev.mat',...
%     'run_4_rev.mat'};
% data_fname = {'run_9_rev.mat',...
%     'run_10_rev.mat'};
data_fname = {'run_1_rev.mat',...
    'run_4_rev.mat',...
    'run_5_rev.mat',...
    'run_6_rev.mat',...
    'run_9_rev.mat',...
    'run_10_rev.mat'};
path_dir = 'C:\Users\David\data\bell_momentumspin\pairsource\90\rev\';
path_main = 'C:\Users\David\Dropbox\PhD\diagrams\data_analysis\k_space\full_analysis';

b_flip_mJ = true;       % flips mJ indices to mJ = 0,1 (default 1,0,-1)


%% main
cd(path_main);
for ii = 1:numel(data_fname)
    data_path = fullfile(path_dir,data_fname{ii});
    
    close all;
    
    % run g2 analysis
    run('collision_atom_distribution.m');    
end
close all;


%% end of script.
t_elapsed = toc;