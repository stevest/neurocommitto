[idx_str,~,~,~,pref_str]=apclusterK(mergedCN_str, ClNo,0)
            [targ_str,~,labels_str] = unique(idx_str);
            histc(labels_str,1:ClNo)'
%%

close all;
rurange = 1;
stop = 5000;
currentCluster = 2;
% load from cluster:
pathprefix = 'X:\Documents\Glia\';
% pathprefix = 'Z:\data\GliaBackup\';
% spathprefix = 'C:\Users\marouli\Documents\';

run.nruns=36;
% RUNS_str = cell(1,run.NC_str);
Sid=1;
fprintf('Loading runs...');
PCcells_str = cell(run.nPC,run.nruns);
% sPCcells_str = cell(run.nPC,run.nruns);
PVcells_str = cell(run.nPC,run.nruns);

ridx = zeros(run.NC_str(Sid),run.nPC,run.nruns);
basalSegments = 5;
dendseg = {};
for stc=currentCluster%:run.NC_str(Sid)
    RUNS_str{1,stc} = cell(run.nPC,run.nruns);
    sRUNS_str{1,stc} = cell(run.nPC,run.nruns);
    for ru = rurange%run.nruns
        ru
        tic;
        S = cell(1,run.nPC);
        for pc=1:run.nPC
            pc
            if (run.ISBINARY)
                try
%     PCcells_str{pc,ru} = ncell(load(sprintf('%s%s/STR_SN%d_ST%d/%d_%d_%d.txt',pathprefix,run.path,run.sn,run.state,stc-1,pc-1,ru-1)),10);
                    PCcells_str{pc,ru} = ncell(nrn_vread(sprintf('%s%s/RND_SN%d_ST%d/%d_%d_%d.bin',pathprefix,run.path,run.sn,run.state,stc-1,pc-1,ru-1),'n'),10);
%                     sPCcells_str{pc,ru} = ncell(load(sprintf('%s%s/STR_SN%d_ST%d/%d_%d_%d.txt',spathprefix,run.path,run.sn,run.state,stc-1,pc-1,ru-1)),10);
                catch err
                    warning('on','verbose')
                    warning(err.message)
                    ridx(stc,pc,ru) = 1; % true, if file not found
                    continue;
                end
            else
            end
            
            
            PCcells_str{pc,ru}.clusterID = run.labels_str(pc,Sid);
%             sPCcells_str{pc,ru}.clusterID = run.labels_str(pc,Sid);
%             [S{pc},~,~] = findUPstates(PCcells_str{pc,ru}.mv(run.stimend*run.dt:run.dt:end),4, 10, -66, 3000 );
            
            PCcells_str{pc,ru}.spikes = single(PCcells_str{pc,ru}.spikes);
%             sPCcells_str{pc,ru}.spikes = single(sPCcells_str{pc,ru}.spikes);
            if (sum(PCcells_str{pc,ru}.spikes>run.stimend)/((stop-run.stimend)/1000)) > 7
                PCcells_str{pc,ru}.persistentActivity = 1;
            else
                PCcells_str{pc,ru}.persistentActivity = 0;
            end
%             PCcells_str{pc,ru}.mv = [];
%             if ~isempty(S{pc});
%                 PCcells_str{pc,ru}.persistentActivity = 1;
% %                 sPCcells_str{pc,ru}.persistentActivity = 1;
%             else
%                 PCcells_str{pc,ru}.persistentActivity = 0;
% %                 sPCcells_str{pc,ru}.persistentActivity = 0;
%             end
%             PCcells_str{pc,ru} = PCcells_str{pc,ru} ;
%             sPCcells_str{pc,ru} = sPCcells_str{pc,ru} ;
            % Load delays PV2PC ( for each target)
%             delayPV2PC_str{pc,ru} = load(sprintf('%s%s/STR_SN%d_ST%d_gaba/delayPV2PC_trg_%d_runs_%d.txt',pathprefix,run.path,run.sn,run.state,pc-1,ru-1));
            
        end
        for pv = 76:100%76:88
            pv
            PVcells_str{pv,ru} = ncell(nrn_vread(sprintf('%s%s/RND_SN%d_ST%d/%d_%d_%d.bin',pathprefix,run.path,run.sn,run.state,stc-1,pv-1,ru-1),'n'),10);
        end
        runtime = toc - tic
    end
    RUNS_str{1,stc} = PCcells_str(:,:);
%     sRUNS_str{1,stc} = sPCcells_str(:,:);

end

% Save memory!
clear PCcells_str;
% clear sPCcells_str;

fprintf('DONE!\n');

steps = run.dt;
bin = 100;
start = run.stimstart;

bins = start:bin:stop;
freqBinsCell = cell(1,36);
CV = zeros(run.nPC,run.NC_str,36);
% Fano = zeros(run.nPC,run.NC_str,36);
ISI = cell(run.nPC,run.NC_str,36);
st = start+1:bin:stop-bin+1;
sp = start+bin:bin:stop;
% spikeBins = zeros(run.nPC,5000,36);
spikeBins = struct();
spikeBinsCluster = struct();
spikeBinsNoCluster = struct();
vCell = cell(1,36);
% GABAaCell = cell(1,36);
% GABAbCell = cell(1,36);
% icaCell = cell(1,36);
% caiCell = cell(1,36);
% inmdaCell = cell(1,36);
clstidx = [];
for k = 1:run.NC_str
    clstidx = [clstidx ; find(run.labels_rnd==k)];
end

if (length(rurange)>1)
    for k=1:run.nPC
        spikeBins(k).spikes = logical(zeros(7,stop));
    end
    clusterCells = find(run.labels_str==currentCluster)';
    noClusterCells = find(run.labels_str~=currentCluster)';
    for k=1:length(clusterCells)
        spikeBinsCluster(k).spikes = logical(zeros(7,stop));
    end

    for k=1:length(noClusterCells)
        spikeBinsNoCluster(k).spikes = logical(zeros(7,stop));
    end
end

% rurange = [2,4:9]

for ru = rurange
    freqBins = zeros(run.nPC,length(bins),run.NC_str);
    vBins = zeros(run.nPC,length(st),run.NC_str);
    for stc = currentCluster%:run.NC_str
        for k = 1:length(clstidx)
            
            spikes = RUNS_str{1,stc}{clstidx(k),ru}.spikes;
            if (length(rurange)>1)
%                 spikeBins(clstidx(k),round(spikes),ru) = 1;
%                 spikes(spikes<start) = [];
                if any(ismember(clusterCells,clstidx(k)))
                    clOvrh = cumsum(run.cellsPerCluster_str(1:currentCluster-1));
                    spikeBinsCluster(k-(clOvrh(end))).spikes(ru,round(spikes)) = 1;
                else

                    spikeBinsNoCluster(k).spikes(ru,round(spikes)) = 1;
                end
                spikeBins(clstidx(k)).spikes(ru,round(spikes)) = 1;
            end
            
            freqhisto = (histc(spikes, bins) / (bin/1000));
            if ~isempty(spikes)
                CV(clstidx(k),stc,ru) = std(diff(spikes)) / mean(diff(spikes));
                ISI{clstidx(k),stc,ru} = diff(spikes(spikes>1500));
            end
            if ~isempty(freqhisto)
                freqBins(clstidx(k),:,stc) = freqhisto ;
            end
            v = RUNS_str{1,stc}{clstidx(k),ru}.mv;
            v = v(1:steps:end);
            for kk = 1:length(st)
                vBins(clstidx(k),kk,stc) = mean(v(st(kk):sp(kk)));
            end
        end
    end
%     sum(freqBins(:))
    freqBinsCell{ru} = freqBins;
    vCell{ru} = vBins;
end


cl=stc
win = bin;
wst = start+1:win:stop-win+1;
wsp = start+win:win:stop;
cAllFF = cell(5,4);
cAllv = cell(5,4);
for ru = rurange
    j = floor((ru-1)/6)+1;
    k = mod(ru-1,6)+1;
    cAllFF{j,k} = freqBinsCell{ru}(clstidx,:,cl);
    cAllv{j,k} = vCell{ru}(clstidx,:,cl);
end
AllFF = cell2mat(cAllFF);
Allv = cell2mat(cAllv);

figure();imagesc(AllFF);hold on;
title('Firing Frequency');
set(gca,'YDir','normal');
set(gca,'Xticklabel',run.stimend+bin*5 : bin*5 :stop-1);
colorbar


cells = find(run.labels_str==stc)'
cellslen = length(cells);
figure;hold on;
cm = hsv(cellslen);
tmp2 = [];
tmp = cell(1,cellslen);
histo = cell(1,cellslen);


for k=1:cellslen
    k
    for l = 1:length(ru)
        tmp{k} = [tmp{k}, ISI{cells(k),stc,ru}];
    end
%     tmp2 = [tmp2, tmp];
    histo{k} = histc(tmp{k}, 0:1:60);

%     maxVal(max(histo{k})>maxVal) = max(histo{k});
%     histoconv = conv(histo{k},fspecial('gaussian',[1 10], 2),'same');
%  histo{k} = (histoconv./max(histoconv))*max(histo{k}) ;
     plot(histo{k});
%     plot(conv(histo{k},fspecial('gaussian',[1 10], 2),'same'));
end
% maxVal = max(cell2mat(cellfun(@max, histo, 'Uniformoutput',false)))

figure;hold on;
for c=1:75 %find(run.labels_str==stc)'
 plot(RUNS_str{1,stc}{c,ru}.mv);hold on;
%  pause;cla;
end
figure;hold on;
for c=76:100
plot(PVcells_str{c,ru}.mv(1:10:end));hold on;
% pause;cla;
end
figure;hold on;
for c=1:75
% 1000/mean(diff(PVcells_str{c,ru}.spikes))
plot(c,length(RUNS_str{1,stc}{c,ru}.spikes) / (stop/1000),'*b')
end
% figure;hold on;
for c=76:100
% 1000/mean(diff(PVcells_str{c,ru}.spikes))
plot(c-75,length(PVcells_str{c,ru}.spikes) / (stop/1000),'*r')
end
