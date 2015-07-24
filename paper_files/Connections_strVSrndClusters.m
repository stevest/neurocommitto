
%Stimulate more cells each time to check network response:
clear all;clc;
pathprefix = 'Z:/data/GliaBackup/BigNetwork_750/';
nStimCells = 50;
run.cellsPerCluster_str
% run.cellsPerCluster_rnd

% Re-run clustering if necessary.
[CN_str] = m_commonNeighbors(run.state_str(1:nPC,1:nPC));
mergedCN_str = CN_str .* run.state_str(1:nPC,1:nPC);
labels_str={};
prefvals = -100:10:30
for pref = 1:length(prefvals)
        NC_str=[];
        [idx_str,~,~,~]=apcluster(mergedCN_str,prefvals(pref),'dampfact',0.9, ...
        'convits',200,'maxits',2000,'nonoise');
        [targ_str,~,labels_str{pref}] = unique(idx_str);
        [within_str(pref),between_str(pref),WB_str(pref)] = calculateWithinBetween(run.state_str(1:nPC,1:nPC),labels_str{pref});
end
plot([within_str;between_str]')

[CN_rnd] = m_commonNeighbors(run.state_rnd(1:nPC,1:nPC));
mergedCN_rnd = CN_rnd .* run.state_rnd(1:nPC,1:nPC);
labels_rnd={};
for pref = 1:length(prefvals)
[idx_rnd,~,~,~]=apcluster(mergedCN_rnd,prefvals(pref),'dampfact',0.9, ...
        'convits',200,'maxits',2000,'nonoise');
        [targ_rnd,~,labels_rnd{pref}] = unique(idx_rnd);
        [within_rnd(pref),between_rnd(pref),WB_rnd(pref)] = calculateWithinBetween(run.state_rnd(1:nPC,1:nPC),labels_rnd{pref});
end
plot([within_rnd;between_rnd]')


