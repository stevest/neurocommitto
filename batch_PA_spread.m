
for aa = 17
    close all; clearvars -except aa; 
    randn('seed',aa);
%     addpath(genpath('~/Documents/MATLAB/'));
%     addpath(genpath('~/Documents/GitHub/prefrontal-micro/'));
%     addpath(genpath('~/Documents/GitHub/neurocommitto/'));
    mypath = '/srv/userdata/HOMES/stefanos/Documents/prefrontal-micro/experiment/network/';
    mypath2 = '/srv/userdata/HOMES/stefanos/Documents/prefrontal-micro/';
    
    exprun=10;
    nPC = 75; %216;
%     clstr = 1; % Cluster type
    
    connBinsPC2PC = 20:30:500;
    connProbsPC2PC = 0.25 .* exp(-0.006 * connBinsPC2PC);
    connBinsPV2PC = 0:20:500;
    connProbsPV2PC = linspace(1,0,length(connBinsPV2PC)) * 0.7;
    connBinsCB2PC = 0:20:500;
    connProbsCB2PC = normpdf(connBinsCB2PC,150,90) *10;
    
    recipBinsPC2PC = 20:30:500;
    recipProbsPC2PC = 0.12 .* exp(-0.006 * recipBinsPC2PC);
    % Update probabilities from Perin et al. (Figure S4)
    NconnBins = [0,1,2,3];
    NincomingProbs = [0.1, 0.2, 0.25, 0.4] ;
    NoutgoingProbs = [0.12, 0.25, 0.24, 0.2] ;
    
    %  ---- Initialize Pyramidal cells ----
    PCsomata = CreateRandomNetwork(nPC, 200, 3);
    distPC2PC = generateDistanceMat(PCsomata', 0);
    
    % Pyramidals connect to all
    PC2PC = zeros(nPC,nPC,2);
    PC2PC(:,:,1) = create_graph_CN(distPC2PC,connBinsPC2PC, connProbsPC2PC, ...%0.07  %0.2 dinei PA
        recipBinsPC2PC, recipProbsPC2PC,NconnBins,NincomingProbs,NoutgoingProbs);
    
    prob_conn_rnd_joint = sum(sum(tril(PC2PC(:,:,1) | PC2PC(:,:,1)') )) / ((numel(PC2PC(:,:,1))-nPC)/2); % individual probability
    prob_conn_rnd_int = sum(sum(PC2PC(:,:,1))) / (numel(PC2PC(:,:,1))-nPC); % joint probability
    
    PC2PC(:,:,2) = create_graph_WS(distPC2PC,prob_conn_rnd_joint,'joint'); % 1.0=Random graph, 0.0=Watts-Strogatz graph
    
    
    [mcccs,~,cccs] = clust_coeff( PC2PC(:,:,1) );
    [mcccr,~,cccr] = clust_coeff( PC2PC(:,:,2) );
    [mcccr,mcccs]
    
    sum(sum(PC2PC(:,:,1) & PC2PC(:,:,1)'))
    sum(sum(PC2PC(:,:,2) & PC2PC(:,:,2)'))

    % routines to calculate and plot the degrees and clustering
    % distribution:
    %     c_m = hot(6);
    %     for i=1:5
    %         p_h(i) = plot(1:5:100,histc(degrees(PC2PC(:,:,i)),1:5:100),'color',c_m(i,:),'linewidth',3);hold on;
    %
    %         set(get(get(p_h(i),'Annotation'),'LegendInformation'),...
    %     'IconDisplayStyle','on');
    %     end
    %     legend('Perin','Distance','Clustered','intermediate','Random');
    
    
    % ---- affinity propagation similarity measure: ----
    
    % Find nearest neighbors of Pyramidals only
    [CNi_str,CNo_str] = m_commonNeighbors(PC2PC(:,:,1));
    [CNi_rnd,CNo_rnd] = m_commonNeighbors(PC2PC(:,:,2));
    mergedCN_str = [CNi_str + CNo_str] .* PC2PC(:,:,1);
    mergedCN_rnd = [CNi_rnd + CNo_rnd] .* PC2PC(:,:,2);
    mergedCN_strT = mergedCN_str';
    mergedCN_str(logical(tril(ones(size(mergedCN_str)),-1))) = mergedCN_strT(logical(tril(ones(size(mergedCN_str)),-1)));
    mergedCN_rndT = mergedCN_rnd';
    mergedCN_rnd(logical(tril(ones(size(mergedCN_rnd)),-1))) = mergedCN_strT(logical(tril(ones(size(mergedCN_rnd)),-1)));
    
    %perform affinity propagation algorithm:
    Sid=1;
    [idx_str,~,~,~]=apclustermex(mergedCN_str,median(mergedCN_str,2)*0.5 );
    [idx_rnd,~,~,~]=apclustermex(mergedCN_rnd,median(mergedCN_rnd,2)*0.5 );
    NC_str(Sid) = length(unique(idx_str));
    [~,~,labels_str(:,Sid)] = unique(idx_str);
    NC_rnd(Sid) = length(unique(idx_rnd));
    [~,~,labels_rnd(:,Sid)] = unique(idx_rnd);
    %     [NC,labels,Sid] = runAffinityPropagation(onlyCN);
    
    
    %how many cells in each cluster?
    cellsPerCluster_str = {histc(labels_str(:,Sid),1:NC_str(Sid))'};
    cellsPerCluster_rnd = {histc(labels_rnd(:,Sid),1:NC_rnd(Sid))'};
    figure('name','Structured')
    bar(cellsPerCluster_str{:})
    figure('name','Random')
    bar(cellsPerCluster_rnd{:})
    
    %  ---- Rest of the cells ----
    nPV = round(nPC*13/75);
    nCB = round(nPC*6/75);
    nCR = round(nPC*6/75);
    if(nPV==0)% NEURON ISSUE:
        nPV = 1;
    end
    if(nCB == 0)
        nCB = 1;
    end
    if(nCR == 0)
        nCR=1;
    end
    nAllCells = sum([nPC,nPV,nCB,nCR]);
    AllCells = [nPC , nPV , nCB , nCR, nAllCells];
    
    PVsomata = CreateCubeNetworkPV(250, nPV); % 226.6 per mm squared (?!?) paper?
    CBsomata = CreateCubeNetworkPV(200, nCB);
    CRsomata = CreateCubeNetworkPV(200, nCR);
    
    % Distance of each interneuron from Pyramidal and connect based on
    % probability from Yuste 2011:
    distPV2PC = distancePV2PC(PVsomata,PCsomata);
    distCB2PC = distancePV2PC(CBsomata,PCsomata);
    
    ConnMatPV2PC = connectPV2PC(distPV2PC,connBinsPV2PC,connProbsPV2PC);
    ConnMatCB2PC = connectCB2PC(distCB2PC,connBinsCB2PC,connProbsCB2PC);
    
    % Connection of FS interneurons to PCs as in :
    % Cortical inhibitory cell types differentially form intralaminar
    % and interlaminar subnetworks with excitatory neurons, Otsuka Takeshi
    % Kawaguchi Yasuo, 2009
    ConnMatPC2PV = connectPC2PV(distPV2PC')  ;
    % Test if the PC2PV reciprocals are 20% of PC2PV pairs:
    ctr = 10000;
    while ctr && (1.5 < abs(((sum(sum((ConnMatPC2PV & ConnMatPV2PC')))*100) / numel(ConnMatPC2PV)) - 20))
        % if reciprocals are more/less, move connections to reach 20%:
        foundIt = 0;
        while ~foundIt
            ip = ceil(rand(1)*size(ConnMatPC2PV,1));
            jp = ceil(rand(1)*size(ConnMatPC2PV,2));
            if ConnMatPC2PV(ip,jp) && ~ConnMatPV2PC(jp,ip)
                foundIt = 1;
            end
        end
        foundIt = 0;
        while ~foundIt
            ii = ceil(rand(1)*size(ConnMatPC2PV,1));
            ji = ceil(rand(1)*size(ConnMatPC2PV,2));
            if ~ConnMatPC2PV(ii,ji) && ConnMatPV2PC(ji,ii)
                foundIt = 1;
            end
        end
        % move the projection of the PC to form reciprocal with the PV:
        ConnMatPC2PV(ip,jp) = 0;
        ConnMatPC2PV(ii,ji) = 1;
        ctr = ctr - 1 ;
    end
    % Set indices for ease of mind:
    pc = size(PCsomata,1);
    pv = size(PCsomata,1) + size(PVsomata,1);
    cb = size(PCsomata,1) + size(PVsomata,1) + size(CBsomata,1);
    cr = size(PCsomata,1) + size(PVsomata,1) + size(CBsomata,1) + size(CRsomata,1);
    
    % Populate the final all-to-all connectivity matrix:
    AllConnMat = zeros(nAllCells);
        
    %         % Pyramidals to all types of interneurons:
    %         AllConnMat(1:pc,pc+1:end) = 1; % Connected to all
    % PCs to PVs
    AllConnMat(1:pc,pc+1:pv) = ConnMatPC2PV;
    AllConnMat(1:pc,pv+1:end) = 1; % Connected to all interneurons
    % PVs connect to all other PVs + autapses (need to gap junctions?)
    AllConnMat(pc+1:pv,pc+1:pv) = 1;
    % PVs to PCs based on above connectivity (Yuste 2011):
    AllConnMat(pc+1:pv,1:pc) = ConnMatPV2PC;
    % CB only connect to PC (Xenia) na to psa3w...
    AllConnMat(pv+1:cb,1:pc) = ConnMatCB2PC;
    % CR only connect to PC and CB (Xenia) na to psa3w...
    AllConnMat(cb+1:cr,1:pc) = 1;
    AllConnMat(cb+1:cr,pv+1:cb) = 1;
    
    AllWeightMat = ones(size(AllConnMat));
    
    
    % Generate stimulation data:
    nIncoming=[130,130];
    randShift = 20;
    nStim = 20;
    tstop = 1000;
    spikesDend = {};
    spikesApic = {};
    for c=1:pc
        i=1;
        while i <= nIncoming(1)
            spikesDend(c,i)={linspace(0,tstop,nStim)};
            spikesDend{c,i} = round(spikesDend{c,i} + ((rand(1,nStim)-0.5)*randShift));
            spikesDend{c,i} = spikesDend{c,i}(spikesDend{c,i}>0);
            spikesDend{c,i} = spikesDend{c,i}(spikesDend{c,i}<tstop);
            spikesDend{c,i} = sort((spikesDend{c,i}),'ascend');
            if(length(unique(spikesDend{c,i})) ~= length(spikesDend{c,i}) )
                continue;
            end
            i = i +1;
        end
        i=1;
        while i <= nIncoming(2)
            spikesApic(c,i)={linspace(0,tstop,nStim)};
            spikesApic{c,i} = round(spikesApic{c,i} + ((rand(1,nStim)-0.5)*randShift));
            spikesApic{c,i} = spikesApic{c,i}(spikesApic{c,i}>0);
            spikesApic{c,i} = spikesApic{c,i}(spikesApic{c,i}<tstop);
            spikesApic{c,i} = sort((spikesApic{c,i}),'ascend');
            if(length(unique(spikesApic{c,i})) ~= length(spikesApic{c,i}) )
                continue;
            end
            i = i +1;
        end
    end
    
%     % GABAa inhibition in pyramidals' somata:
%     % approximate No of PV in range of each pyramidal connected
%     % propabilisticaly with distance as in R.Yuste, 2011
%     PVdummies = rand(nPV*5*5*3,3);  % 226.6 per mm squared (?!?) paper?
%     PVdummies = PVdummies .* repmat([750,750,450],length(PVdummies),1);
%     distPVd2PC = distancePV2PC(PVdummies,[375,375,225]);
%     distPVd2PC(distPVd2PC >= 400) = 0;
%     ConnMatPVd2PC = connectPV2PC(distPVd2PC,connBinsPV2PC,connProbsPV2PC);
%     sum(ConnMatPVd2PC) % No of inhibitory GABAa synapses to a pyramidal

    % Generate background inhibition data:
    nIncoming=[396];% hardcoded to avoid run to run variability
    randShift = 1000/6;
    nStim = 30; % ~6Hz as in Alain Destexhe & Denis Pare, 1999
    tstop = 5000;
    spikesGABAa = {};
    for c=1:pc
        i=1;
        while i <= nIncoming(1)
            spikesGABAa(c,i)={linspace(0,tstop,nStim)};
            spikesGABAa{c,i} = round(spikesGABAa{c,i} + ((rand(1,nStim)-0.5)*randShift));
            spikesGABAa{c,i} = spikesGABAa{c,i}(spikesGABAa{c,i}>0);
            spikesGABAa{c,i} = spikesGABAa{c,i}(spikesGABAa{c,i}<tstop);
            spikesGABAa{c,i} = sort((spikesGABAa{c,i}),'ascend');
            if(length(unique(spikesGABAa{c,i})) ~= length(spikesGABAa{c,i}) )
                continue;
            end
            i = i +1;
        end
    end
    
    % run for the structured network and for the random network as
    % well:
    
    RUNS_str = {};
    RUNS_rnd = {};
    
%     export network parameters for that experiment:
    cd(mypath);
    system(sprintf('mkdir experiment_%d',aa));
    exportNetworkPositions(PCsomata,PVsomata,sprintf('%sexperiment_%d/',mypath,aa));
    % Export pyramidal's clusters ID:
    exportNetworkCluster([idx_str,idx_rnd],sprintf('%sexperiment_%d/',mypath,aa));
    % Export stimulation for each pyramidal:
    exportNetworkStimulation(spikesDend,spikesApic,sprintf('%sexperiment_%d/',mypath,aa));
    exportNetworkInhibition(spikesGABAa,sprintf('%sexperiment_%d/',mypath,aa));
    
    save(sprintf('%sexperiment_%d/experiment.mat',mypath,aa),'-v7.3');
    
    for t=1:NC_str(Sid)
        fprintf('@experiment %d, cluster %d... ',aa, t);
        % Isolate one cluster:
        [r,~,~] = find(labels_str(:,Sid) == t);
        StimVect_str = zeros(nPC,1);
        StimVect_str(r) = 1;   
%         StimVect_str = ones(nPC,1); %stimulate WHOLE NETWORK
        
        AllConnMat_str = AllConnMat;
        AllConnMat_str(1:pc,1:pc) = PC2PC(:,:,1);
        AllConnMat_str(sub2ind([nPC,nPC],[1:pc],[1:pc])) = 1;   
        
        AllWeightMat_str = AllWeightMat;
        
        % Kayaguchi MONO GIA TA RECIPROCAL!
        % OXI TA STO IDIO CLUSTER
        % Telika mallon afta pou einai sto idio cluster 
        % H prepei na koita3w kai tis dyo periptwseis...
        for ec = 1:NC_str(Sid)
            clid = find(labels_str(:,Sid) == ec);
            AllWeightMat_str(clid,clid) = 1.3;
        end
        
        AllWeightMat_str = AllWeightMat_str .* AllConnMat_str ;
        
        %  ---- Export parameters to NEURON ----
        exportNetworkParameters(AllCells,AllConnMat_str,AllWeightMat_str,sprintf('%sexperiment_%d/',mypath,aa));
        exportStimulationParameters(AllCells,StimVect_str,sprintf('%sexperiment_%d/',mypath,aa));
        
        fprintf('Running NEURON... ');
        unix(sprintf('kill -9 `ps -ef | grep nrniv | grep -v grep | awk ''{print$2}''`'));
        tic;
%         [nrnStatus,nrnCmdOut{t}] = unix('./parallel_pfc');
        EXPERIMENT = 1;  % 0 = random, 1= structured gia to $#$#% NEURON
%         [nrnStatus,nrnCmdOut{t}] = unix(sprintf('mpiexec -np 8 %smechanism_complex/x86_64/special -mpi -nobanner -c "PARALLEL=1" -c "SIMPLIFIED=0" -c "CLUSTER_ID=%d" -c "EXPERIMENT=%d"  -c "AA=%d" final.hoc',mypath2,t,EXPERIMENT,aa));
        [nrnStatus,nrnCmdOut{t}] = unix(sprintf('mpirun -n 24 -mca btl ^openib %smechanism_complex/x86_64/special -mpi -nobanner -c "PARALLEL=1" -c "SIMPLIFIED=0" -c "CLUSTER_ID=%d" -c "EXPERIMENT=%d"  -c "AA=%d" final.hoc',mypath2,t,EXPERIMENT,aa));
        runtime = toc
        
        if( ~findstr('Success',nrnCmdOut{t}))
            continue;
        end
        
        
        fprintf('DONE!\n');
    end % for t
    
    tt = (cellsPerCluster_str{:}); %% unique
    for t=1:length(tt) %1:NC_rnd(Sid)
        fprintf('@experiment %d, cluster %d... ',aa, t);
        % Isolate one cluster:
%         [r,~,~] = find(labels_rnd(:,Sid) == t);
        StimVect_rnd = zeros(nPC,1);
%         StimVect_rnd(r) = 1;   
%         StimVect_rnd = ones(nPC,1); %stimulate WHOLE NETWORK
%         sum(labels_str == 2)
%         StimVect_rnd = 
        
        AllConnMat_rnd = AllConnMat;
        AllConnMat_rnd(1:pc,1:pc) = PC2PC(:,:,2);
        AllConnMat_rnd(sub2ind([nPC,nPC],[1:pc],[1:pc])) = 1;
        
        AllWeightMat_rnd = AllWeightMat;
        
        [~,isort] = sort(degrees(AllConnMat_rnd(1:pc,1:pc)),'descend')
        StimVect_rnd(isort(1:tt(t))) = 1;
        find(StimVect_rnd)
        
        % Kayaguchi MONO GIA TA RECIPROCAL!
        % OXI TA STO IDIO CLUSTER
        % Morissima , 2011 , journal of neuroscience, ta weights!
        %         for ec = 1:NC(Sid)
        %             clid = find(labels(:,Sid) == ec);
        %             AllWeightMat(clid,clid) = 1.2;
        %         end

        %  ---- Export parameters to NEURON ----

        exportNetworkParameters(AllCells,AllConnMat_rnd,AllWeightMat_rnd,sprintf('%sexperiment_%d/',mypath,aa));
        exportStimulationParameters(AllCells,StimVect_rnd,sprintf('%sexperiment_%d/',mypath,aa));
%         exportNetworkStimulation(AllCells,NetStimVect_rnd,sprintf('%sexperiment_%d/',mypath,aa))
        
        fprintf('Running NEURON... ');
        unix(sprintf('kill -9 `ps -ef | grep nrniv | grep -v grep | awk ''{print$2}''`'));
        tic;
%         [nrnStatus,nrnCmdOut{t}] = unix('./parallel_pfc');
        EXPERIMENT = 0;  % 0 = random, 1= structured gia to $#$#% NEURON
        [nrnStatus,nrnCmdOut{t}] = unix(sprintf('mpirun -n 24 -mca btl ^openib %smechanism_complex/x86_64/special -mpi -nobanner -c "PARALLEL=1" -c "SIMPLIFIED=0" -c "CLUSTER_ID=%d" -c "EXPERIMENT=%d" -c "AA=%d" final.hoc',mypath2,t,EXPERIMENT,aa));
        runtime = toc
        
        if( ~findstr('Success',nrnCmdOut{t}))
            continue;
        end
        
        
        fprintf('DONE!\n');
        
        %         else
        %             continue;
        %         end
    end % for t
    
    
    
end

NC_rnd=NC_str;

%% Load cluster info:
mypath = '/srv/userdata/HOMES/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/';
mypath2 = '/srv/userdata/HOMES/stefanos/Documents/GitHub/prefrontal-micro/';
idx_str = load(sprintf('%sexperiment_%d/networkPyramidalClustersStructured.txt',mypath,aa));
idx_rnd = load(sprintf('%sexperiment_%d/networkPyramidalClustersRandom.txt',mypath,aa));
[~,~,labels_str] = unique(idx_str);
[~,~,labels_rnd] = unique(idx_rnd);
Sid = 1;
exprun = 10;
pc = 75;
NC_str = length(unique(labels_str));
NC_rnd = length(unique(labels_rnd));

PCsomata = load(sprintf('%sexperiment_%d/networkPositionsPyramidals.txt',mypath,aa));

%%
fprintf('Loading runs...');
PCcells_str = {};
for t=1:NC_str(Sid)
    for ru = 1:exprun
        for c=1:pc
            mycell = ncell(load(sprintf('%sexperiment_%d/%s/%d_%d_%d.txt',mypath,aa,'STR',t,c-1,ru-1)),10);
            mycell.clusterID = labels_str(c,Sid);
            mycell.position = PCsomata(c,1:3);
            PCcells_str{c,ru}=mycell.hasPersistent(1000,8,mycell.tstop-1001); % paper?
        end
    end
    RUNS_str{1,t} = PCcells_str(:,:);
end
fprintf('DONE!\n');

fprintf('Loading runs...');
PCcells_rnd = {};
for t=1:NC_rnd(Sid)
    for ru = 1:exprun
        for c=1:pc
            mycell = ncell(load(sprintf('%sexperiment_%d/%s/%d_%d_%d.txt',mypath,aa,'RND',t,c-1,ru-1)),10);
            mycell.clusterID = labels_rnd(c,Sid);
            mycell.position = PCsomata(c,1:3);
            PCcells_rnd{c,ru}=mycell.hasPersistent(1000,8,mycell.tstop-1001); % paper?
        end
    end
    RUNS_rnd{1,t} = PCcells_rnd(:,:);
end
fprintf('DONE!\n');
    
%     save(sprintf('../../../%d_N%d_R%d_C%d.mat',aa,clstr,exprun,labels_str(Sid)),'-v7.3');
%% Plot Connectivity matrices
figure;hold on;
[row,col] = find(AllConnMat_str(1:pc,1:pc))
scatter(row,col,'o', 'fill');hold on;



%%
for i=1:NC(Sid)
    clustID(:,i) = (labels(:,Sid) == i);
    CPA{i} = cellfun(@(x) x.persistentActivity,RUNS{1,i}(clustID(:,i),:));
    NCPA{i} = cellfun(@(x) x.persistentActivity,RUNS{1,i}(~clustID(:,i),:));
end

% generate color maps:
cm = rand(NC(Sid),3);
idcs= mat2cell(1:NC(Sid),[1],ones(NC(Sid),1));

figure();hold on;
cellfun(@(x,i) plot(histc(sum(x,1),0:1:nPC),'color',cm(i,:)), CPA,idcs,'uniformoutput',false);

%% Plot-a-cluster
close all;

clu = 1 ;% what cluster is stimulated
ru = 1;% what run of this cluster

CellsPerClust=zeros(NC_str(Sid),1);
PAPerClust=zeros(NC_str(Sid),1);


for cl=1:NC_str(Sid)
    figure(); hold on;
    set(gcf,'numbertitle','off','name', sprintf('Cluster #%d',cl) );
    for c=1:pc
        if(RUNS_str{1,clu}{c,ru}.clusterID == cl)
            %                         figure;
            plot( RUNS_str{1,clu}{c,ru}.mv,'linewidth',1,'color',rand(1,3));
            CellsPerClust(cl) = CellsPerClust(cl)+1;
            if(RUNS_str{1,clu}{c,ru}.persistentActivity)
                PAPerClust(cl) = PAPerClust(cl) +1;
            end
        end
    end
    axis([0,length(RUNS_str{1,clu}{c,ru}.mv),-70,70 ]);hold on;
end

% bar graph
figure();
bar(CellsPerClust,'b');hold on;
bar(PAPerClust,'r');

% Extract some stats:

%% Plot Random clusters
close all;

clu = 2 ;% what cluster is stimulated
ru = 2;% what run of this cluster

CellsPerClust=zeros(NC_str(Sid),1);
PAPerClust=zeros(NC_str(Sid),1);


for cl=1:NC_str(Sid)
    figure(); hold on;
    set(gcf,'numbertitle','off','name', sprintf('Cluster #%d',cl) );
    for c=1:pc
        if(RUNS_str{1,clu}{c,ru}.clusterID == cl)
            %                         figure;
            plot( RUNS_str{1,clu}{c,ru}.mv,'linewidth',1,'color',rand(1,3));
            CellsPerClust(cl) = CellsPerClust(cl)+1;
            if(RUNS_str{1,clu}{c,ru}.persistentActivity)
                PAPerClust(cl) = PAPerClust(cl) +1;
            end
        end
    end
    axis([0,length(RUNS_str{1,clu}{c,ru}.mv),-70,70 ]);hold on;
end

% bar graph
figure();
bar(CellsPerClust,'b');hold on;
bar(PAPerClust,'r');

% Extract some stats:
%% Probability of recruitment:
close all;
PAP_mean = [];
PAP_isfinite_sum=[];
PAP_recruit_mean = [];

for dataset=1:2
% load dataset:
if dataset == 1
    NC = NC_str(Sid);
    RUNS = RUNS_str;
    cpc = histc(labels_str(:,Sid),1:NC_str(Sid)) ;
else 
    NC = NC_rnd(Sid);
    RUNS = RUNS_rnd;
    cpc = histc(labels_rnd(:,Sid),1:NC_rnd(Sid)) ;
end


for ru = 1:10 % what run of this cluster
    
    CellsPerClust=zeros(NC(Sid));
    PAPerClust=zeros(NC(Sid));
    
    for clu =1:NC(Sid)
        for cl=1:NC(Sid)
            for c=1:pc
                if(RUNS{1,clu}{c,ru}.clusterID == cl)
                    if(RUNS{1,clu}{c,ru}.persistentActivity)
                        PAPerClust(cl,clu) = PAPerClust(cl,clu) +1;
                    end
                end
            end
        end
    end
    PAP(:,:,ru) = (PAPerClust.*(~eye(length(PAPerClust))));
end

% % Figure #Active/#inCluster:
% PAP(:,:,1) ./ repmat(cpc,1,NC(Sid))
% 
% % Figure #Absolut cells:
% PAP(:,:,1)
% 
% % Figure #Clusters Recruited:
% sum(PAP>=1,3)/10


PAPSD = reshape(sum(PAP),[NC(Sid),10,1])';
mean(PAPSD)
std(PAPSD)
cpc = histc(labels_str(:,Sid),1:NC_str(Sid)) ;
PAP_recruit = [];
if dataset == 1
    PAP_mean(:,:,1) = mean(PAP,3);
    PAP_isfinite_sum(:,:,1) = sum(PAP>0,3);
    for iii = 1:exprun
       PAP_recruit(:,:,iii) = PAP(:,:,iii)   ./ repmat(cpc,1,NC_str) ;
    end
    PAP_recruit_mean(:,:,1) = mean(PAP_recruit,3);
else 
    PAP_mean(:,:,2) = mean(PAP,3);
    PAP_isfinite_sum(:,:,2) = sum(PAP>0,3);
    for iii = 1:exprun
       PAP_recruit(:,:,iii) = PAP(:,:,iii)   ./ repmat(cpc,1,NC_str) ;
    end
    PAP_recruit_mean(:,:,2) = mean(PAP_recruit,3);
end

end %for datasets

% ccm = jet( ceil(max(max(max(PCFF_mean)))*100)  ) ; 
% % ccm = ccm(1:100:end,:);
% figure('Name','Structured Firing Frequency' );
% imagesc(PCFF_mean(:,:,1));
% colormap(ccm(1:round(max(max(PCFF_mean(:,:,1)))*100),:));
% colorbar;
% set(gcf,'Color',[1,1,1]);

ccm = jet( ceil(max(max(max(PAP_recruit_mean)))*100)  ) ; 
figure('Name','Recruitement Structured' );
imagesc(PAP_recruit_mean(:,:,1));
colormap(ccm(1:round(max(max(PAP_recruit_mean(:,:,1)))*100),:));
colorbar;
set(gcf,'Color',[1,1,1]);

figure('Name','Recruitement Random' );
imagesc(PAP_recruit_mean(:,:,2));
colormap(ccm(1:round(max(max(PAP_recruit_mean(:,:,2)))*100),:));
colorbar;
set(gcf,'Color',[1,1,1]);


% figure('Name','Recruitement Structured' );
% imagesc(PAP_mean(:,:,1)/10);
% colorbar;
% 
% figure('Name','Recruitement Random' );
% imagesc(PAP_mean(:,:,2)/10);
% colorbar;


figure('Name','Histo Recruitement' );
ph = bar( 0:0.02:0.2,histc(reshape(PAP_recruit_mean(:,:,1),1,[]),0:0.02:0.2 )' );hold on;
ch = get(ph,'Children');
set(ch,'FaceAlpha',0.6, 'FaceColor', 'r') ;
ph = bar( 0:0.02:0.2,histc(reshape(PAP_recruit_mean(:,:,2),1,[]),0:0.02:0.2 )' );hold on;
ch = get(ph,'Children');
set(ch,'FaceAlpha',0.6, 'FaceColor', 'b') ;
set(gcf,'Color',[1,1,1]);

%% Panos simple:

close all;

for dataset=1:2
% load dataset:
if dataset == 1
    NC = NC_str(Sid);
    RUNS = RUNS_str;
else 
    NC = NC_rnd(Sid);
    RUNS = RUNS_rnd;
end

PS={};
for ru = 1:exprun % what run of this cluster
    ResponcePerClust=cell(NC(Sid));
    
    for stimclu =1:NC(Sid)
        for cl=1:NC(Sid)
            for c=1:pc
                if(RUNS{1,stimclu}{c,ru}.clusterID == cl)
                        responce = RUNS{1,stimclu}{c,ru}.mv(10001:end);
                        ResponcePerClust{cl,stimclu}(c,1) = {find( ([0;diff(sign(diff(responce)))<0;0] & [sign(responce)==1]) )./10 };
                end
            end
            % Panos simple:
            PS{cl,stimclu,ru} =  ResponcePerClust{cl,stimclu}(~cellfun(@isempty,ResponcePerClust{cl,stimclu} ) ) ;
        end
    end
end

% define windowing function:
inWin = @(x,a,b)(a<=x) & (x<=b)


% PER CLUSTER (between clusters):
windows=[];
windows(1,:) = 0:200:4800;
windows(2,:) = 200:200:5000;
WFFPC = {};
WFFIC = {};
for i=1:length(windows)
    minWin = windows(1,i);
    maxWin = windows(2,i);
    windowing_PC = @(a)sum(cellfun(@sum ,cellfun(@(x) inWin(x,minWin,maxWin), a,'uniformoutput',false)))  ;
    windowing_IC = @(a)cellfun(@sum ,cellfun(@(x) inWin(x,minWin,maxWin), a,'uniformoutput',false)) ;
    WFFPC{1,i} = cellfun(windowing_PC, PS);
    WFFIC{1,i} = cellfun(windowing_IC, PS,'uniformoutput',false);
end



%% PLOTING:

% Plot per cluster:
figure; hold on;
run = 3
blah=[];
stimCluster = 3 ;
for select = 1:NC(Sid) %select cluster
    blah(select,:) = cellfun(@(x) x(select,stimCluster,run) ,WFFPC);
end
normfact = mean(mean(blah));
blah = blah - normfact;
set(gca, 'ColorOrder', rand(NC(Sid),3), 'NextPlot', 'replacechildren');
plot(blah'); hold on;
plot(sum(blah,1),'Color', 'k' , 'linewidth',2) ;


%% Plot inside cluster:
figure; hold on;
run = 3  ;
theCluster = 1 ;
stimCluster = 1 ;
blah=[];
blah= cell2mat(cellfun(@(x) x(theCluster,stimCluster,run) ,WFFIC));
normfact = mean(mean(blah));
blah = blah - normfact;
set(gca, 'ColorOrder', rand(size(blah,1),3), 'NextPlot', 'replacechildren');
plot(blah'); hold on;
plot(sum(blah,1),'Color', 'k' ) ;


%% Plot individual cells:
figure; hold on;
run = 3  ;
stimCluster = 1 ;
blah=[];
for theCluster=1:NC(Sid)
    blah = [blah ; cell2mat(cellfun(@(x) x(theCluster,stimCluster,run) ,WFFIC)) ] ;
end
normfact = mean(mean(blah));
blah = blah - normfact;
set(gca, 'ColorOrder', rand(size(blah,1),3), 'NextPlot', 'replacechildren');
plot(blah'); hold on;
plot(sum(blah,1),'Color', 'k' ) ;




%% Firing Frequency:
close all;
PCFF_mean = [];

for dataset=1:2
% load dataset:
if dataset == 1
    NC = NC_str(Sid);
    RUNS = RUNS_str;
else 
    NC = NC_rnd(Sid);
    RUNS = RUNS_rnd;
end

FF = [];
TMP=[];
for ru = 1:exprun % what run of this cluster
    FFPerClust=cell(NC(Sid));
    PAPerClust=zeros(NC(Sid));
    
    for stimclu =1:NC(Sid)
        for cl=1:NC(Sid)
            for c=1:pc
                if(RUNS{1,stimclu}{c,ru}.clusterID == cl)
                    if(RUNS{1,stimclu}{c,ru}.persistentActivity)
                        %                         PAPerClust(cl,stimclu) = PAPerClust(cl,stimclu) +1;
                        responce = RUNS{1,stimclu}{c,ru}.mv(10001:end);
                        FFPerClust{cl,stimclu} = [FFPerClust{cl,stimclu} ,length(find( ([0;diff(sign(diff(responce)))<0;0] & [sign(responce)==1]) ))];
                    end
                end
            end
        end
    end
    %     PAP(:,:,ru) = (PAPerClust.*(~eye(length(PAPerClust))));
    %     FFPC(:,:,ru) = cellfun(@(x,i) i*(nanmean(x)),FFPerClust,mat2cell(cellfun(@(x) ~isempty(x), FFPerClust),ones(NC(Sid),1),ones(NC(Sid),1)) );
    
    
    TMP = cellfun(@(x) nanmean(x),FFPerClust );
    TMP(~isfinite(TMP))  = 0;
    FF(:,:,ru) = TMP./ ((RUNS{1,stimclu}{c,ru}.tstop-1001)/1000)  ; % Get Freq, no num of spikes
end

FFnonstim = FF;

FFstim = [];
for i=1:exprun
    FFstim = [FFstim;diag(FFnonstim(:,:,i)) ];
    FFnonstim(:,:,i) = FFnonstim(:,:,i).*(~eye(NC(Sid)));
end
FFnonstim =FFnonstim(find(FFnonstim))  ;

% non stimulated
mean(FFnonstim)
std(FFnonstim)

% in stimulated cluster:
mean(FFstim)
std(FFstim)

% Per Cluster:
if dataset == 1
    PCFF_mean(:,:,1) = mean(FF,3);
else 
    PCFF_mean(:,:,2) = mean(FF,3);
end


% ttest between the same cluster:
for i=1:NC(Sid)
    for j=1:NC(Sid)
        [h(i,j),p(i,j)] = ttest2(reshape(FF(i,:,:),1,[]),reshape(FF(j,:,:),1,[]))
    end
end

end% for datasets

ccm = jet( ceil(max(max(max(PCFF_mean)))*100)  ) ; 
% ccm = ccm(1:100:end,:);
figure('Name','Structured Firing Frequency' );
imagesc(PCFF_mean(:,:,1));
colormap(ccm(1:round(max(max(PCFF_mean(:,:,1)))*100),:));
colorbar;
set(gcf,'Color',[1,1,1]);

figure('Name','Random Firing Frequency' );
imagesc(PCFF_mean(:,:,2));
colormap(ccm(1:round(max(max(PCFF_mean(:,:,2)))*100),:));
colorbar;
set(gcf,'Color',[1,1,1]);

figure('Name','Histo Firing Frequency' );
ph = bar( 0:0.5:ceil(max(max(max(PCFF_mean)))), histc(reshape(PCFF_mean(:,:,1),1,[]),0:0.5:ceil(max(max(max(PCFF_mean)))) )' );hold on;
ch = get(ph,'Children');
set(ch,'FaceAlpha',0.6, 'FaceColor', 'b') ;
ph = bar( 0:0.5:ceil(max(max(max(PCFF_mean)))), histc(reshape(PCFF_mean(:,:,2),1,[]),0:0.5:ceil(max(max(max(PCFF_mean)))) )' );hold on;
ch = get(ph,'Children');
set(ch,'FaceAlpha',0.6, 'FaceColor', 'r') ;
set(gcf,'Color',[1,1,1]);

figure;
bar(cpc);
set(gcf,'Color',[1,1,1]);


% for i=1:NC(Sid)
% %     PCFFstim(1:exprun) = reshape(FF(i,i,:),1,[])
% %     PCFFnonstim(1:NC(Sid),1:exprun) = reshape(FF(:,i,:),NC(Sid),[]);
% %     PCFFnonstim(i,:) = [];
% %     
% %     
% %     
% %     %Nonstimulated:
% %     mean(PCFFnonstim,2)
% end

%% Coefficient of Variation:

close all;
ISI_mean = [];

for dataset=1:2
% load dataset:
if dataset == 1
    NC = NC_str(Sid);
    RUNS = RUNS_str;
else 
    NC = NC_rnd(Sid);
    RUNS = RUNS_rnd;
end

ISI = [];
for ru = 1:exprun % what run of this cluster
    ISIPerClust=cell(NC(Sid));
%     PAPerClust=zeros(NC(Sid));
    
    for stimclu =1:NC(Sid)
        for cl=1:NC(Sid)
            for c=1:pc
                if(RUNS{1,stimclu}{c,ru}.clusterID == cl)
                    if(RUNS{1,stimclu}{c,ru}.persistentActivity)
                        whichones = RUNS{1,stimclu}{c,ru}.spikes >= 1001; % do not count the stimulus spiking
                        ISIPerClust{cl,stimclu,c} = [diff(RUNS{1,stimclu}{c,ru}.spikes(whichones))'];
                    end
                end
            end
        end
    end
    
    for cl=1:NC(Sid)
        for stimclu =1:NC(Sid)
            TMP = reshape(ISIPerClust(cl,stimclu,:),[],1,1);
            TMP = TMP(~cellfun(@isempty,TMP)); % each cell seperately (per cell)
            ISI(cl,stimclu,ru) = nanmean(cellfun(@std ,TMP) ./ cellfun(@mean ,TMP)); % mean of individual CV s (per cluster)
        end
    end
%     ISI(:,:,ru) = ISIPerClust;%./4  ;
end

if dataset == 1
    ISI_mean(:,:,1) = (nanmean(ISI,3));
else 
    ISI_mean(:,:,2) = (nanmean(ISI,3));
end


end % for dataset

ccm = jet( ceil(max(max(max(ISI_mean)))*100)  ) ; 
% ccm = ccm(1:100:end,:);
figure('Name','Structured CV' );
imagesc(ISI_mean(:,:,1));
colormap(ccm(1:round(max(max(ISI_mean(:,:,1)))*100),:));
colorbar;
set(gcf,'Color',[1,1,1]);

figure('Name','Random CV' );
imagesc(ISI_mean(:,:,2));
colormap(ccm(1:round(max(max(ISI_mean(:,:,2)))*100),:));
colorbar;
set(gcf,'Color',[1,1,1]);

figure('Name','Histo CV' );
ph = bar( 0:0.05:ceil(max(max(max(ISI_mean)))), histc(reshape(ISI_mean(:,:,1),1,[]),0:0.05:ceil(max(max(max(ISI_mean)))) )' );hold on;
ch = get(ph,'Children');
set(ch,'FaceAlpha',0.6, 'FaceColor', 'b') ;
ph = bar( 0:0.05:ceil(max(max(max(ISI_mean)))), histc(reshape(ISI_mean(:,:,2),1,[]),0:0.05:ceil(max(max(max(ISI_mean)))) )' );hold on;
ch = get(ph,'Children');
set(ch,'FaceAlpha',0.6, 'FaceColor', 'r') ;
set(gcf,'Color',[1,1,1]);

% % Per cluster:
% for i=1:NC(Sid)
%     for j=1:NC(Sid)
%         for k = 1:exprun
%             CV(i,j,k) = std(ISI{i,j,k}) / mean(ISI{i,j,k});
%         end
%     end
% end
% 
% mean(CV,3)
% 
% 
% % for non stimulated condition (per cluster):
% ISIPCtmp = ISI;
% 
% for k = 1:exprun
%     for i=1:NC(Sid)
%         ISIPCtmp{i,i,k} = [];
%     end
% end
% 
% ISIinclust={};
% for i=1:NC(Sid)
%     ISIinclust{i}=[];
%     for k=1:exprun
%         ISIinclust{i} = [ISIinclust{i},ISIPCtmp{i,:,k}];
%     end
% end
% figure('Name','Non stimulated');
% bar(cellfun(@std,ISIinclust)./cellfun(@mean,ISIinclust))
% 
% % for the stimulated condition: (per cluster)
% ISIPCdiag = {};
% for i=1:NC(Sid)
%     ISIPCdiag{i} = [];
%     for k = 1:exprun
%         ISIPCdiag{i} = [ISIPCdiag{i},ISI{i,i,k}] ;
%     end
% end
% figure('Name','Stimulated');
% bar(cellfun(@std,ISIPCdiag)./cellfun(@mean,ISIPCdiag))

%% Van Rossum / Spiky multi neuron distance:
close all;

para.tmin = 0;
para.tmax = 4000;
para.dts =  10;
para.dtm = 10;
para.select_measures=[0 1     0 0     0 0     0 0];            % Select order of measures

PairwiseSynch_mean = [];

for dataset=1:2
% load dataset:
if dataset == 1
    NC = NC_str(Sid);
    RUNS = RUNS_str;
else 
    NC = NC_rnd(Sid);
    RUNS = RUNS_rnd;
end


SynchPerClust = cell(NC(Sid),NC(Sid),exprun,pc);
PairwiseSynch = [];
for ru = 1:exprun % what run of this cluster
    for stimclu =1:NC(Sid)
        for cl=1:NC(Sid)
            for c=1:pc
                if(RUNS{1,stimclu}{c,ru}.clusterID == cl)
                    if(RUNS{1,stimclu}{c,ru}.persistentActivity)
                        %                         PAPerClust(cl,stimclu) = PAPerClust(cl,stimclu) +1;
                        whichones = RUNS{1,stimclu}{c,ru}.spikes >= 1001;
                        SynchPerClust{cl,stimclu,ru,c} =RUNS{1,stimclu}{c,ru}.spikes(whichones);
                    end
                end
            end
        end
    end
    
%     for ru=1:exprun
        for stimclu =1:NC(Sid)
            STIM = reshape(SynchPerClust(stimclu,stimclu,ru,:),[],1,1);
            STIM = STIM(~cellfun(@isempty,STIM));
            if size(STIM,1) == 0
                continue;
            end
            for cl=1:NC(Sid)
                if(cl ~= stimclu)
                    min_elements=[];
                    tmpels=[];
                    min_rndperm_els=[];
                NONSTIM = reshape(SynchPerClust(cl,stimclu,ru,:),[],1,1);
                NONSTIM = NONSTIM(~cellfun(@isempty,NONSTIM));
                if size(NONSTIM,1) == 0
                    continue;
                end
%                 min_elements = min(size(STIM,1), size(NONSTIM,1));
%                 tmpels = randperm(size(STIM,1));
%                 min_rndperm_els(1,:) = tmpels(1:min_elements);
%                 tmpels = randperm(size(NONSTIM,1));
%                 min_rndperm_els(2,:) = tmpels(1:min_elements);
%                 PairwiseSynch(cl,stimclu,ru) = MNvanRossum(STIM(min_rndperm_els(1,:)),NONSTIM(min_rndperm_els(2,:)),20,1)  ;
                % take all the pairwise comparisons (so, I can use
                % different synchronicity measure from Spiky to cross
                % correlation:
                cumres = 0;
                for comb_stim=1:size(STIM,1)
                    for comb_nost=1:size(NONSTIM,1)
                        results = SPIKY_no_plot_f_distances_MEX({STIM{comb_stim}';NONSTIM{comb_nost}'},para);
                        cumres = cumres + results.overall_dissimilarities ;
                    end
                end
                PairwiseSynch(cl,stimclu,ru) = cumres / (size(STIM,1)*size(NONSTIM,1));
                
                end
                
            end
        end
%     end
end
if dataset == 1
    PairwiseSynch_mean(:,:,1) = nanmean(PairwiseSynch,3);
else 
    PairwiseSynch_mean(:,:,2) = nanmean(PairwiseSynch,3);
end


end % for dataset

ccm = jet( ceil(max(max(max(PairwiseSynch_mean)))*1000)  ) ; 
% ccm = ccm(1:100:end,:);
figure('Name','Structured Synchronicity' );
imagesc(PairwiseSynch_mean(:,:,1));
colormap(ccm(1:round(max(max(PairwiseSynch_mean(:,:,1)))*1000),:));
colorbar;
set(gcf,'Color',[1,1,1]);

figure('Name','Random Synchronicity' );
imagesc(PairwiseSynch_mean(:,:,2));
colormap(ccm(1:round(max(max(PairwiseSynch_mean(:,:,2)))*1000),:));
colorbar;
set(gcf,'Color',[1,1,1]);

figure('Name','Histo Synchronicity' );
ph = bar( 0:0.005:(max(max(max(PairwiseSynch_mean)))), histc(reshape(PairwiseSynch_mean(:,:,1),1,[]),0:0.005:(max(max(max(PairwiseSynch_mean)))) )' );hold on;
ch = get(ph,'Children');
set(ch,'FaceAlpha',0.6, 'FaceColor', 'b') ;
ph = bar( 0:0.005:(max(max(max(PairwiseSynch_mean)))), histc(reshape(PairwiseSynch_mean(:,:,2),1,[]),0:0.005:(max(max(max(PairwiseSynch_mean)))) )' );hold on;
ch = get(ph,'Children');
set(ch,'FaceAlpha',0.6, 'FaceColor', 'r') ;
set(gcf,'Color',[1,1,1]);

%% Cross correlation:
close all;
PerClust = cell(NC(Sid),NC(Sid),exprun,nPC);
for ru = 1:exprun % what run of this cluster
    for stimclu =1:NC(Sid)
        for cl=1:NC(Sid)
            for c=1:nPC
                if(RUNS{1,stimclu}{c,ru}.clusterID == cl)
                    if(RUNS{1,stimclu}{c,ru}.persistentActivity)
                        %                         PAPerClust(cl,stimclu) = PAPerClust(cl,stimclu) +1;
                        tmpSpikes = RUNS{1,stimclu}{c,ru}.spikes;
                        PerClust{cl,stimclu,ru,c} = tmpSpikes(tmpSpikes > 1001);
                    end
                end
            end
        end
    end
end

% Init SPIKY:
para.tmin = 0;
para.tmax = 5000;
para.dts =  10;
para.dtm = 10;
para.select_measures=[0 1     0 0     0 0     0 0];            % Select order of measures
% para.select_measures=[1];            % Select order of measures


% Check synchronicity of cells belonging to the stimulated cluster:
synchInsideClust={};
synchInsideClustclean={};
overal_dissim_stimulated={};
for i=1:NC(Sid)
    synchInsideClust{i}=[];
    for k=1:exprun
        synchInsideClust(i,k) = {reshape(PerClust(i,i,k,reshape(~cellfun(@isempty,PerClust(i,i,k,:)),1,[])),1,[])};
    end
    %concatinate all the cells of this cluster/run to one pool:
    % get idx of runs containing more than 1 cell firing after stimuly
    % (PA):
    goodIdx = cellfun(@(x) x>1,cellfun(@(x) length(x),synchInsideClust(i,:),'uniformoutput',false));
    synchInsideClustclean(i,1:sum(goodIdx)) = synchInsideClust(i,goodIdx) ;
    
    % feed SPIKY:
    for k=1:sum(goodIdx)
        results = SPIKY_no_plot_f_distances_MEX(synchInsideClustclean{i,k},para);
        overal_dissim_stimulated{i,k} = results.overall_dissimilarities ;
    end
    
    %     synchInsideClust{i,~cellfun(@isempty,synchInsideClust(i,:))}
end
% mean synchronicity of all the stimulated clusters:
mean_synch_stim = mean(cell2mat(overal_dissim_stimulated(cellfun(@(x) ~isempty(x), overal_dissim_stimulated))));
std_synch_stim = std(cell2mat(overal_dissim_stimulated(cellfun(@(x) ~isempty(x), overal_dissim_stimulated))));

% Check synchronicity of cells belonging to the non-stimulated cluster:
synchOutsideClust={};
syncnOutsideClustclean={};
overal_dissim_nostim={};
for i=1:NC(Sid)
    synchOutsideClust{i}=[];
    ctrk=1;
    for j=1:NC(Sid)
        for k=1:exprun
            if i ~= j
                synchOutsideClust(i,ctrk) = {reshape(PerClust(i,j,k,reshape(~cellfun(@isempty,PerClust(i,j,k,:)),1,[])),1,[])};
                ctrk = ctrk+1;
            end
        end
    end
    %concatinate all the cells of this cluster/run to one pool:
    % get idx of runs containing more than 1 cell firing after stimuly
    % (PA):
    goodIdx = cellfun(@(x) x>1,cellfun(@(x) length(x),synchOutsideClust(i,:),'uniformoutput',false));
    syncnOutsideClustclean(i,1:sum(goodIdx)) = synchOutsideClust(i,goodIdx) ;
    
    % feed SPIKY:
    for k=1:sum(goodIdx)
        results = SPIKY_no_plot_f_distances_MEX(syncnOutsideClustclean{i,k},para);
        overal_dissim_nostim{i,k} = results.overall_dissimilarities ;
    end
    
    %     synchInsideClust{i,~cellfun(@isempty,synchInsideClust(i,:))}
end
% mean synchronicity of all the stimulated clusters:
mean_synch_nostim = mean(cell2mat(overal_dissim_nostim(cellfun(@(x) ~isempty(x), overal_dissim_nostim))));
std_synch_nostim = std(cell2mat(overal_dissim_nostim(cellfun(@(x) ~isempty(x), overal_dissim_nostim))));


% Check synchronicity between stimulated and non stimulated clusters
% (pairwise):
% trains_a = {};
%     trains_b ={};
% for i=1:11
%     trains_a{i} = ceil(rand(ceil(rand(1)*20),1)*100) ;
%     trains_b{i} = ceil(rand(ceil(rand(1)*20),1)*100) ;
% end
% trains_a = trains_a';
% trains_b = trains_b';
% output = MNvanRossum(trains_a,trains_b,0.9,0.1)

syncnPairInClust={};
syncnPairInClustClean={};
syncnPairOutClust={};
syncnPairOutClustClean={};
overal_dissim_pairwise={};
for i=1:NC(Sid)
    for j=1:NC(Sid)
        % get spikes from stimulated cluster:
        if(i ~= j)
            %         fprintf('taking from ii %d\n',i)
            syncnPairInClust{j}=[];
            ctrk=1;
            for k=1:exprun
                syncnPairInClust(j,ctrk) = {reshape(PerClust(i,i,k,reshape(~cellfun(@isempty,PerClust(i,i,k,:)),1,[])),1,[])} ;
                ctrk = ctrk+1;
            end
            goodIdxIn = cellfun(@(x) x>1,cellfun(@(x) length(x),syncnPairInClust(j,:),'uniformoutput',false));
            syncnPairInClustClean(i,1:sum(goodIdxIn)) = syncnPairInClust(j,goodIdxIn) ;
            
            
            
            % get spikes from all other clusters in different cells:
            syncnPairOutClust{j}=[];
            ctrk=1;
            for k=1:exprun
                if i ~= j
                    syncnPairOutClust(j,ctrk) = {reshape(PerClust(i,j,k,reshape(~cellfun(@isempty,PerClust(i,j,k,:)),1,[])),1,[])};
                    ctrk = ctrk+1;
                end
            end
            
            goodIdxOut = cellfun(@(x) x>1,cellfun(@(x) length(x),syncnPairOutClust(j,:),'uniformoutput',false));
            syncnPairOutClustClean(i,1:sum(goodIdxOut)) = syncnPairOutClust(j,goodIdxOut) ;
            
            min_elements = [];
            min_rndperm_els=[];
            if( ~isempty(syncnPairOutClustClean) && ~isempty(syncnPairInClustClean) )
                min_elements = min(length(syncnPairInClustClean{1}), length(syncnPairOutClustClean{1}));
                tmpels = randperm(length(syncnPairInClustClean{1}));
                min_rndperm_els(1,:) = tmpels(1:min_elements);
                tmpels = randperm(length(syncnPairOutClustClean{1}));
                min_rndperm_els(2,:) = tmpels(1:min_elements);
                
                min_elements = min(length(syncnPairInClustClean{1}), length(syncnPairOutClustClean{1}));
                randperm(min_elements);
                overal_dissim_pairwise{i,j} = MNvanRossum(syncnPairInClustClean{1}(min_rndperm_els(1,:))',syncnPairOutClustClean{1}(min_rndperm_els(2,:) )',0.9,0.1)  ;
                %             fprintf('comparing with jj %d\n',j)
            end
            
        end
    end
end

overal_dissim_pairwise(cellfun(@isempty,overal_dissim_pairwise)) = {0};
overal_dissim_pairwise = cell2mat(overal_dissim_pairwise);
%%
scatter3(PCsomata(:,1),PCsomata(:,2), PCsomata(:,3), 20, [0,0,1], 'fill', 'o');
axis equal; hold on;