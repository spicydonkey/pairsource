% Monte-carlo simulation of multimode thermal statistics
% DKS 2020

%% config
% thermal distribution
p_thermal = @(mu,x) (mu.^x)./((mu+1).^(x+1));

Nsamp = 1e1;            % number of samples
nShot = 1e4;            % number of shots

M = 1e3;                % number of modes
x1_mean = 3e-3;         % mean atom number in single mode

x1_supp = 0:10;         % support for single mode (cutoff)

X_ed = (0:(M*x1_mean*10))-0.5;
X_ct = edge2cent(X_ed);



%% Monte-carlo simulation
p1 = thermalpdf(x1_supp,x1_mean);        % single mode prob vector


progressbar('MC sim');

p_X = NaN(Nsamp,length(X_ct));
for ii = 1:Nsamp
    x1s = reshape(randsample(x1_supp,nShot*M,true,p1),[nShot,M]);
    tX = sum(x1s,2);            % this detected number X
    
    % get prob dist from this sample
    tp_X = nhist(tX,{X_ed});
    tp_X = tp_X/sum(tp_X);      % normalise probability
    
    p_X(ii,:) = tp_X;
    
    progressbar(ii/Nsamp);
end


%% statistic
p_X_mean = mean(p_X,1);
p_X_std = std(p_X,0,1);


%% plot
h_mc = figure('name','MC_sim');
hold on;

ploterr(X_ct,p_X_mean,[],p_X_std,'hhxy',0);

% poisson
Y_poiss = poisspdf(X_ct,x1_mean*M);
plot(X_ct,Y_poiss,'k--');

ax = gca;
ax.YScale = 'log';

