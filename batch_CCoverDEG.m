close all; clear all; clc;
mypath = '/srv/userdata/HOMES/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/';
cd(mypath);

exprun=1:1;
nPC = 64;%216;

for t=exprun
    fprintf('@experiment %d... ',t);
    %  ---- Initialize Pyramidal cells ----
    % Initialize somata positions
    PCsomata = CreateRandomNetwork(nPC, 200, 3);
    distPC2PC = generateDistanceMat(PCsomata', 0);
    
    % Pyramidals connect to all
    checkCC = -1;
    while checkCC<=0 || checkCC>1 || sum(sum(PC2PC))<3
        PC2PC = create_graph_WS(distPC2PC,rand(1),0.99); % 1.0=Random graph, 0.0=Watts-Strogatz graph
        checkCC = clust_coeff(PC2PC);
    end
    Clustering(t) = checkCC;
    Degree(t) = mean(degrees(PC2PC));
    %  ---- Rest of the cells ----
    nPV = round(nPC*17/65);
    nCB = round(nPC*9/65);
    nCR = round(nPC*9/65);
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
    
    % Populate the final all-to-all connectivity matrix:
    AllConnMat = zeros(nAllCells);
    % Set indies for ease of mind:
    pc = size(PCsomata,1);
    pv = size(PCsomata,1) + size(PVsomata,1);
    cb = size(PCsomata,1) + size(PVsomata,1) + size(CBsomata,1);
    cr = size(PCsomata,1) + size(PVsomata,1) + size(CBsomata,1) + size(CRsomata,1);
    
    PC2PC(sub2ind([nPC,nPC],[1:pc],[1:pc])) = 1; % autapses in Pyramidals
    AllConnMat(1:pc,1:pc) = PC2PC;
    
    % Pyramidals to all types of interneurons:
    AllConnMat(1:pc,pc+1:end) = 1; % Connected to all
    % PVs connect to all other PVs + autapses (need to gap junctions?)
    AllConnMat(pc+1:pv,pc+1:pv) = 1;
    % PVs to PCs based on above connectivity (Yuste 2011):
    AllConnMat(pc+1:pv,1:pc) = 1;
    % CB only connect to PC (Xenia) na to psa3w...
    AllConnMat(pv+1:cb,1:pc) = 1;
    % CR only connect to PC and CB (Xenia) na to psa3w...
    AllConnMat(cb+1:cr,1:pc) = 1;
    AllConnMat(cb+1:cr,pv+1:cb) = 1;
    
    %  ---- Export parameters to NEURON ----
    cd(mypath);
    unix('rm -rf multi_core');
    unix('mkdir multi_core');
    exportNetworkParameters(AllCells,AllConnMat,mypath);
    exportStimulationParameters(AllCells,ones(nPC,1),mypath);
    
    fprintf('Running NEURON... ');
    [nrnStatus,nrnCmdOut{t}] = unix('./parallel_pfc');
    if(nrnStatus~=0)
        continue;
    end
    
    fprintf('Loading runs... ');
    for ru = 1:5
        for c=1:nPC
            mycell = ncell(load(sprintf('multi_core/%d_%d.txt',c-1,ru-1)),10);
            PCcells{c,ru}=mycell.hasPersistent(1000,9,4000); % paper?
        end
    end
    RUNS{t} = PCcells(:,:);
    Percentages(t) = sum(all(cellfun(@(x) x.persistentActivity,RUNS{1,t}))) * 100 / size(RUNS{1,t},2);
    fprintf('Saving results... ');
    save('batch_CC_DEG.mat');
    fprintf('DONE!\n');
end % for t

%%
for ii=1:5
    plot(RUNS{1}{1,ii}.mV)
    pause;
    cla
end

%%

% for k=1:size(RUNS,2)
%         % count if ALL cells of the cluster have persistent:
%         Percentages(k) = sum(all(cellfun(@(x) x.persistentActivity,RUNS{1,k}))) * 100 / size(RUNS{1,k},2);
% 
% end
