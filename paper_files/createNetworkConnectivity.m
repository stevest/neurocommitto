% Stefanos, 2014 
% Built a structured and random network.  
close all; clear all; clc

% Random seed
aa=1;
rng('default') 
rng(aa)

%Set paths
% mypath = '/home/cluster/papoutsi/Desktop/STEFANOS/nassi/experiment/network/matlab_files';
% mypath2 = '/home/cluster/papoutsi/Desktop/STEFANOS/nassi/';
% addpath('/home/cluster/papoutsi/Desktop/STEFANOS/nassi/experiment/network/matlab_files/adapt_apV3')
% addpath('/home/cluster/papoutsi/Desktop/STEFANOS/nassi/experiment/network/matlab_files/code')

% Number of neurons
nPC = 75;
nPV = round(nPC*25/75);
if(nPV==0)
    nPV = 1;
end
nAllCells = sum([nPC,nPV]);
AllCells = [nPC , nPV , nAllCells];

%  ---- 3D positions ----
cube_dimensions=200;
% 3-d points in space. Arguments: #of cells, max seperation distance.
[PCsomata, distPC2PC]= CreateRandomNetwork(nPC, cube_dimensions);

% 3d points for PV. Returns 3d points and distances from PCcells  
% Nassi -Max seperation Distance as PC-PC
[PVsomata, distPV2PC] = CreateCubeNetworkPV(cube_dimensions, nPV, PCsomata);

% figure;
% scatter3(PCsomata(:,1),PCsomata(:,2),PCsomata(:,3)); 
% hold on
% scatter3(PVsomata(:,1),PVsomata(:,2),PVsomata(:,3), 'r');

%% Connect PV-PC neurons 
% Distance of each interneuron from Pyramidal and connect based on
% probability from Yuste 2011:    
% gia to connectivity twn PV 2 PC o Packer et al., 2011:
% elegxei gia ena region gyrw apo to ka8e PC. Epomenws ta histograms einai
% swsta gia ta PV2PC (ka8ws briskontai sto idio layer). Oso gia ta
% CB(SOM)2PC pou briskontai se diaforetika layers mporoume na eikasoume
% oti:
% * Apo to sxhma twn SOM (Packer et al., 2013), oso metakineisai sto 
% transverse layer oi tuft dendrites enws PC 8a exoun tis idies pi8anotites
% gia overlap opos exoun kai ta PC pou briskontai sto idio layer. Epomenws
% metraw san distance CB2PC mono to distance sto transverse plane kai
% xrisimopoiw tis pi8anotites pou dinei o Packer et al., 2013, Fig4D
% Packer, A. M., & Yuste, R. (2011). Dense, unspecific connectivity of 
% neocortical parvalbumin-positive interneurons: a canonical microcircuit 
% for inhibition? The Journal of Neuroscience : The Official Journal of the
% Society for Neuroscience, 31(37), 1326071. 
% doi:10.1523/JNEUROSCI.3131-11.2011

% Packer, A., McConnell, D., Fino, E., & Yuste, R. (2013). Axo-dendritic 
% overlap and laminar projection can explain interneuron connectivity to 
% pyramidal cells. Cerebral Cortex, (December), 27902802. 
% doi:10.1093/cercor/bhs210

% Connection of FS interneurons to PCs as in :
% Cortical inhibitory cell types differentially form intralaminar
% and interlaminar subnetworks with excitatory neurons, Otsuka Takeshi
% Kawaguchi Yasuo, 2009
% Random, with probability 23%
  
% PV-PC
% Nassi
connProbsPV2PC= @(x) ( 0.4 .* exp(-0.003 * x));

ConnMatPV2PC = connectPV2PC(distPV2PC,connProbsPV2PC);
ConnMatPC2PV = connectPC2PV(distPV2PC');

% Test if the PC2PV reciprocals are 20% of all pairs. Based on Otsuka, Kaywagushi-Figure2B. Rearrange connection accordingly.:
ctr = 10000;
while ctr && (0.1 < abs(((sum(sum((ConnMatPC2PV & ConnMatPV2PC')))*100) / numel(ConnMatPC2PV)) - 20))
    
    % if reciprocals are more/less, move connections to reach 20%:
    foundIt = 0;
    while ~foundIt
        ip = ceil(rand(1)*size(ConnMatPC2PV,1));
        jp = ceil(rand(1)*size(ConnMatPC2PV,2));
%         If there PC-to-PV connection but not PV to PC
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
% Otsuka: No connection 0.7, reciprocal 0.2, PC-IN 3%, IN-PC 7%
Con_prob(1,1)=(sum(sum((ConnMatPC2PV==0 & ConnMatPV2PC'==0)))*100) / numel(ConnMatPC2PV);
Con_prob(2,1)=(sum(sum((ConnMatPC2PV==1 & ConnMatPV2PC'==1)))*100) / numel(ConnMatPC2PV);
Con_prob(3,1)=(sum(sum((ConnMatPC2PV==1 & ConnMatPV2PC'==0)))*100) / numel(ConnMatPC2PV);
Con_prob(4,1)=(sum(sum((ConnMatPC2PV==0 & ConnMatPV2PC'==1)))*100) / numel(ConnMatPC2PV);

%PV-PV. Connection Probability: O.77 (Gibson, Connors, 1999)
%Eh, o Gibson einai layers 4 and 6, eno o Galarreta, Hestrin, 1999 (same
%issue) einai smoatosensory/visual L5:
% 66% electrical coupling <=80? appart.
% 18% GABAergic connections
% 4.5% Reciprocal GABAergic connections
% ConnMatPV2PV=connectPV2PV(nPV);
[ConnMatPV2PV, gapmat]=connectPV2PV(nPV);
%Label the gap junctions for NEURON (pontiako implementation):
start = 2;
ctr = 1;
for k = 1:nPV
    for j = start:nPV
        if gapmat(k,j)
            gapmat(k,j) = ctr;
            ctr = ctr + 2;
        end
    end
end
gtmp = gapmat';
gapmat(logical(tril(ones(nPV),-1))) = gtmp(logical(tril(ones(nPV),-1)));

% Set indices for ease of mind:
pc = size(PCsomata,1);
pv = size(PCsomata,1) + size(PVsomata,1);
   
% Populate the final all-to-all connectivity matrix. In this matrix rows
% are the sources, columns are the targets.
AllConnMat = zeros(nAllCells);
        
% Pyramidals to interneurons:
% PCs to PVs
AllConnMat(1:pc,pc+1:pv) = ConnMatPC2PV;
% PVs connect to all other PVs + autapses 
AllConnMat(pc+1:pv,pc+1:pv) = connectPV2PV(nPV);
% PVs to PCs based on above connectivity (Yuste 2011):
AllConnMat(pc+1:pv,1:pc) = ConnMatPV2PC;
  
% imagesc(AllConnMat)
  

%% Connect pyramidals
% Distance-dependent connection probabilities.
% PC-PC, Non- reciprocal (Nassi)
% BinsPC2PC = 0:30:400;

% Nassis (fer reference)
% total_prob = @(x) 0.22.*exp(-0.006*x);
% recipProbsPC2PC = @(x) 0.15 .* exp(-0.01* x);
% connProbsPC2PC =@(x) ((0.22.*exp(-0.006*x)) - ( 0.12 .* exp(-0.01* x)))*2;
%%
close all;
total_prob = @(x) 0.22.*exp(-0.0052*x);
recipProbsPC2PC = @(x) 0.12 .* exp(-0.0085* x);
connProbsPC2PC =@(x) ((0.22.*exp(-0.0052*x)) - (0.12 .* exp(-0.0085* x)))*2;


overall = imread('overall.jpg');
reciprocal = imread('reciprocal.jpg');
nonreciprocal = imread('nonreciprocal.jpg');
fact = 665;
figure;imagesc(flipud(overall));hold on;
set(gca,'Ydir','normal');
plot(total_prob(0:250)*fact,'r');

figure;imagesc(flipud(nonreciprocal));hold on;
set(gca,'Ydir','normal');
plot(connProbsPC2PC(0:250)*fact,'k');

figure;imagesc(flipud(reciprocal));hold on;
set(gca,'Ydir','normal');
plot(recipProbsPC2PC(0:250)*fact,'b');
%%
many_nets=1;

PC2PC_str = zeros(nPC,nPC,many_nets);
PC2PC_d   = zeros(nPC,nPC,many_nets);
PC2PC_rnd = zeros(nPC,nPC,many_nets);

for i=1:many_nets
    [PC2PC_d(:,:,i), PC2PC_str(:,:,i), coeffs_global(i, 2:3), coeffs_local(i,2:3), ~, pp(i,2:3)]= create_graph_CN(distPC2PC,connProbsPC2PC,recipProbsPC2PC, total_prob);
    prob_conn_rnd_ind =pp(1,3);%0.13; 
    % na to tre3w polles fores k na kratisw afto pou einai pio konta sto
    % input independent probability...
    [PC2PC_rnd(:,:,i), coeffs_global(i, 1), coeffs_local(i,1), pp(i, 1)] = create_graph_WS(nPC,prob_conn_rnd_ind,'ind'); % 1.0=Random graph, 0.0=Watts-Strogatz graph
end

sum(PC2PC_rnd(:))
sum(PC2PC_d(:))
sum(PC2PC_str(:))
%No of connected pairs is different because in random case the connections
%get distributed so more pairs are connected.
fprintf('Random connected pairs are %.3f%%\n',((sum(sum((PC2PC_rnd | PC2PC_rnd').*triu(ones(nPC),1)))) / ((nPC^2 - nPC)/2)*100));
fprintf('Distance connected pairs are %.3f%%\n',((sum(sum((PC2PC_d | PC2PC_d').*triu(ones(nPC),1)))) / ((nPC^2 - nPC)/2)*100));
fprintf('Structured connected pairs are %.3f%%\n',((sum(sum((PC2PC_str | PC2PC_str').*triu(ones(nPC),1)))) / ((nPC^2 - nPC)/2)*100));

fprintf('Random reciprocal pairs are %.3f%%\n',((sum(sum(PC2PC_rnd & PC2PC_rnd'))/2) / ((nPC^2 - nPC)/2)*100));
fprintf('Distance reciprocal pairs are %.3f%%\n',((sum(sum(PC2PC_d & PC2PC_d'))/2) / ((nPC^2 - nPC)/2)*100));
fprintf('Structured reciprocal pairs are %.3f%%\n',((sum(sum(PC2PC_str & PC2PC_str'))/2) / ((nPC^2 - nPC)/2)*100));

% distPC2PC = states.gparams(1).distPC2PC;
orderedDist = triu(distPC2PC,1);
orderedDist = orderedDist(orderedDist>0);
conn = PC2PC_str;
% conn = states.state_str(1:75,1:75,1);
myrange = 1:10:250;
orderedConn = conn(:,:,end)|conn(:,:,end)';
orderedConn = orderedConn & logical(triu(ones(size(distPC2PC))));
orderedConn = distPC2PC(logical(orderedConn));

distBins = histc( orderedDist, myrange);
bar(distBins)

figure
bar(histc( orderedConn, myrange)./distBins)
set(gca,'XTick',1:length(myrange))
set(gca,'XTickLabel',myrange)


%% Find clusters in networks

for stateID = 1:many_nets
    % Find nearest neighbors.
    [CN_str] = m_commonNeighbors(PC2PC_str(:,:,stateID));
    [ CN_d ] = m_commonNeighbors(PC2PC_d(:,:,stateID));
    [CN_rnd] = m_commonNeighbors(PC2PC_rnd(:,:,stateID));
    
    % keep only the pairs that are connected.
    mergedCN_str = CN_str .* PC2PC_str(:,:,stateID);
    mergedCN_d   =  CN_d  .* PC2PC_d(:,:,stateID);
    mergedCN_rnd = CN_rnd .* PC2PC_rnd(:,:,stateID);
    
    % Affinity propagation algorithm (Frey, Dueck, 2007):
    % Force different # of cluster to get as many clusters as you
    % can from the populations with high min cells per cluster:
    ClNo = 10;
    Flg = 1;
    while Flg
        for i=1:5
            labels_str=[];
            NC_str=[];
            % Nassi, Changed in apclusterK total # of bisections from 20 to 50 (line 75)
            [idx_str,~,~,~,pref_str]=apclusterK(mergedCN_str, ClNo,0)
            [targ_str,~,labels_str] = unique(idx_str);
            if ( min(histc(labels_str,1:ClNo)') > 6)
                Flg = 0;
                ClNo
                break;
            end
        end
        ClNo = ClNo - 1
    end
    Flg = 1;
    while Flg
        for i=1:20
            labels_str=[];
            NC_str=[];
            % Nassi, Changed in apclusterK total # of bisections from 20 to 50 (line 75)
            [idx_d,~,~,~,pref_d]=apclusterK(mergedCN_d, ClNo,0);
            [targ_d,~,labels_d] = unique(idx_d);
            if (min(histc(labels_d,1:ClNo)')> 6)
                Flg = 0;
                ClNo
                break;
            end
        end
        ClNo = ClNo - 1
    end
    Flg = 1;
    while Flg
        for i=1:20
            labels_str=[];
            NC_str=[];
            % Nassi, Changed in apclusterK total # of bisections from 20 to 50 (line 75)
            [idx_rnd,~,~,~,pref_rnd]=apclusterK(mergedCN_rnd, ClNo,0);
            [targ_rnd,~,labels_rnd] = unique(idx_rnd);
            if (min(histc(labels_rnd,1:ClNo)') > 6)
                Flg = 0;
                ClNo
                break;
            end
        end
        ClNo = ClNo - 1
    end
    
    
    cl_labels_str(:,stateID)=labels_str;
    cl_labels_d(:,stateID)=labels_d;
    cl_labels_rnd(:,stateID)=labels_rnd;
    cellsPerCluster_str(stateID) = {histc(labels_str,1:size(targ_str,1))'};
    figure('name','Structured')
    bar(cellsPerCluster_str{:, stateID});
    
    cellsPerCluster_d(stateID) = {histc(labels_d,1:size(targ_d,1))'};
    figure('name','Distance')
    bar(cellsPerCluster_d{:, stateID});
    
    cellsPerCluster_rnd(stateID)= {histc(labels_rnd,1:size(targ_rnd,1))'};
    figure('name','Random')
    bar(cellsPerCluster_rnd{:, stateID});
    
    [within_old,between_old,WB_old] = calculateWithinBetween(states.state_str(1:75,1:75,4),states.gparams(4).labels_str)
    [within_str,between_str,WB_str] = calculateWithinBetween(PC2PC_str,labels_str);
    [within_d,between_d,WB_d] = calculateWithinBetween(PC2PC_d,labels_d);
    [within_rnd,between_rnd,WB_rnd] = calculateWithinBetween(PC2PC_rnd,labels_rnd);
    [within_rnd,within_d,within_str]
    [between_rnd,between_d,between_str]
    
    bar(1,[within_str],'r');hold on;
    bar(2,[between_str],'k');hold on;
    legend({'Intra','Inter'})

    
%     temp_str=[];
%     temp_d=[];
%     temp_rnd=[];
%     for i=1:size(targ_str,1)
%         within=PC2PC_str(find(idx_str==targ_str(i)),find(idx_str==targ_str(i)),stateID);
%         between=PC2PC_str(find(idx_str~=targ_str(i)),find(idx_str~=targ_str(i)),stateID);
%         p_str{stateID,i,1}=sum(sum(within))/((size(within,1)*size(within,1))-size(within,1));
%         p_str{stateID,i,2}=sum(sum( between))/((size( between,1)*size( between,1))-size( between,1));
%         temp_str(i)=p_str{stateID,i,1}- p_str{stateID,i,2};
%     end
%     for i=1:size(targ_d,1)
%         within=PC2PC_d(find(idx_d==targ_d(i)),find(idx_d==targ_d(i)),stateID);
%         between=PC2PC_d(find(idx_d~=targ_d(i)),find(idx_d~=targ_d(i)),stateID);
%         p_d{stateID,i,1}=sum(sum(within))/((size(within,1)*size(within,1))-size(within,1));
%         p_d{stateID,i,2}=sum(sum(between))/((size( between,1)*size( between,1))-size( between,1));
%         temp_d(i)=p_d{stateID,i,1}- p_d{stateID,i,2};
%     end
%     for i=1:size(targ_rnd,1)
%         within=PC2PC_rnd(find(idx_rnd==targ_rnd(i)),find(idx_rnd==targ_rnd(i)),stateID);
%         between=PC2PC_rnd(find(idx_rnd~=targ_rnd(i)),find(idx_rnd~=targ_rnd(i)),stateID);
%         p_rnd{stateID,i,1}=sum(sum(within))/((size(within,1)*size(within,1))-size(within,1));
%         p_rnd{stateID,i,2}=sum(sum( between))/((size( between,1)*size( between,1))-size( between,1));
%         temp_rnd(i)=p_rnd{stateID,i,1}- p_rnd{stateID,i,2};
%     end
    
end
% p_dif(stateID,1)= mean(temp_str);
% p_dif(stateID,2)= mean(temp_d);
% p_dif(stateID,3)= mean(temp_rnd);



%%
% Set indices for ease of mind:
pc = size(PCsomata,1);
pv = size(PCsomata,1) + size(PVsomata,1);
   
% Populate the final all-to-all connectivity matrix. In this matrix rows
% are the sources, columns are the targets.
AllConnMat = zeros(nAllCells);
        
% Pyramidals to interneurons:
% PCs to PVs
AllConnMat(1:pc,pc+1:pv) = ConnMatPC2PV;
% PVs connect to all other PVs + autapses 
AllConnMat(pc+1:pv,pc+1:pv) = connectPV2PV(nPV);
% PVs to PCs based on above connectivity (Yuste 2011):
AllConnMat(pc+1:pv,1:pc) = ConnMatPV2PC;
  
% imagesc(AllConnMat)
  
AllWeightMat = ones(size(AllConnMat));
        
% run for the structured network and for the random network as
% well:

RUNS_str = {};
RUNS_rnd = {};

        
        
        %     export network parameters for that experiment:
        cd(mypath);
        system(sprintf('mkdir experiment_%d',aa));
%         exportNetworkPositions(PCsomata,PVsomata,sprintf('%sexperiment_%d/',mypath,aa));
%         % Export pyramidal's clusters ID:
%         exportNetworkCluster([idx_str,idx_rnd],sprintf('%sexperiment_%d/',mypath,aa));
%         
%         save(sprintf('%sexperiment_%d/experiment.mat',mypath,aa),'-v7.3');

     
        Clusters_str(1:nPC,stateID) = cl_labels_str;

        stimCellsPerCluster = min(min(cellsPerCluster_str{1,stateID}),min(cellsPerCluster_rnd{1,stateID}));
        
    
        for t=1:max(Clusters_str(1:nPC,stateID))
            fprintf('@experiment %d, cluster %d... ',aa, t);
            % Isolate one cluster:
            [r,~,~] = find(Clusters_str(1:nPC,stateID) == t);
            StimVect_str = zeros(nPC,1);
            tmpr = randperm(length(r));
            StimVect_str(r(tmpr(1:stimCellsPerCluster))) = 1;  % 7==NUmber of total cells stimulated in each cluster.
            %         StimVect_str = ones(nPC,1); %stimulate WHOLE NETWORK
            Clusters_stim_str(1:nPC,stateID) = StimVect_str;
            % Autapses:
            AllConnMat_str = AllConnMat;
            AllConnMat_str(1:pc,1:pc) = PC2PC(:,:,1);
            AllConnMat_str(sub2ind([nPC,nPC],[1:pc],[1:pc])) = 1;
            %         AllConnMat_str(1:pc,1:pc) = 1;
            
            AllWeightMat_str = AllWeightMat;
            
            % Kayaguchi : more potent synapses in reciprocal pairs
            % Perin: more potent synapses with more clustering coeff:
            % We try to combine them
            [~,~,cccs] = clust_coeff( AllConnMat_str(1:pc,1:pc) );
            cccs = cccs-min(cccs)
            cccs = cccs/max(cccs)
            weightsexp = exp( [0:0.001:0.3]);
            weightsexp = max(weightsexp) - weightsexp
            weightsexp = weightsexp(end:-1:1);
            plot(weightsexp)
            
            for c = 1:nPC
%                 AllWeightMat_str(1:nPC,c) = 1 + (cccs(c)*1.5);
                    AllWeightMat_str(1:nPC,c) = 1 + (cccs(c)*weightsexp(c)*1.5);
            end
            
            %         AllWeightMat_str = AllWeightMat_str .* AllConnMat_str ;
%             AllWeightMat_str = ones(size(AllConnMat_str));
            
            %  ---- Export parameters to NEURON ----
            exportNetworkParameters(AllCells,AllConnMat_str,AllWeightMat_str,sprintf('%sexperiment_%d/',mypath,aa));
            exportStimulationParameters(AllCells,StimVect_str,sprintf('%sexperiment_%d/',mypath,aa));
            
            fprintf('Running NEURON... ');
            unix(sprintf('kill -9 `ps -ef | grep stefanos | grep nrniv | grep -v grep | awk ''{print$2}''`'));
            tic;
            VCLAMP = 0;
            EXPERIMENT =1;  % 0 = random, 1= structured gia to $#$#% NEURON
            SIMPLIFIED = 1;
            PARALLEL = 1;
            if (SIMPLIFIED)
                neuroncmd = sprintf(['mpirun -n 24 -mca btl ^openib '...
                    '%smechanism_simple/x86_64/special -mpi -nobanner '...
                    '-c "PARALLEL=%d" '...
                    '-c "SIMPLIFIED=%d" '...
                    '-c "CLUSTER_ID=%d" '...
                    '-c "EXPERIMENT=%d" '...
                    '-c "EXPID=%d" '...
                    '-c "AA=%d" '...
                    '-c "VCLAMP=%f" '...
                    'final.hoc'], ...
                    mypath2,PARALLEL,SIMPLIFIED,t,EXPERIMENT,stateID,aa, VCLAMP);
            else
                neuroncmd = sprintf(['mpirun -n 24 -mca btl ^openib '...
                    '%smechanism_complex/x86_64/special -mpi -nobanner '...
                    '-c "PARALLEL=%d" '...
                    '-c "SIMPLIFIED=%d" '...
                    '-c "CLUSTER_ID=%d" '...
                    '-c "EXPERIMENT=%d" '...
                    '-c "EXPID=%d" '...
                    '-c "AA=%d" '...
                    '-c "VCLAMP=%f" '...
                    'final.hoc'], ...
                    mypath2,PARALLEL,SIMPLIFIED,t,EXPERIMENT,stateID,aa, VCLAMP);
            end
            [nrnStatus,nrnCmdOut{t}] = unix(neuroncmd);
            runtime = toc
            
            if( ~findstr('Success',nrnCmdOut{t}))
                continue;
            end
            
            fprintf('DONE!\n');
        end % for t
end
%% Load cluster info:
mypath = '/srv/userdata/HOMES/stefanos/Desktop/prefrontal-micro/experiment/network/';
mypath2 = '~/Desktop/prefrontal-micro/';
load(sprintf('%sexperiment_%d/experiment.mat',mypath,aa));
idx_str = load(sprintf('%sexperiment_%d/networkPyramidalClustersStructured.txt',mypath,aa));
idx_rnd = load(sprintf('%sexperiment_%d/networkPyramidalClustersRandom.txt',mypath,aa));
[~,~,labels_str] = unique(idx_str);
[~,~,labels_rnd] = unique(idx_rnd);
Sid = 1;
exprun = 1;
pc = 75;
stimend = 1500;
NC_str = length(unique(labels_str));
NC_rnd = length(unique(labels_rnd));
PCsomata = load(sprintf('%sexperiment_%d/networkPositionsPyramidals.txt',mypath,aa));

%%
fprintf('Loading runs...');
tic;
%choose structured or random:
% expString = 'RND';
expString = sprintf('STR_%d',stateID);

PCcells_str = {};
PCcells_str_spk = {};
PCcells_str_i = {};
PCcells_str_dist= {};
for t=1:max(Clusters_str(1:nPC,stateID))
    for ru = exprun
        for c=1:pc
            mycell = ncell(load(sprintf('%sexperiment_%d/%s/%d_%d_%d.txt',mypath,aa,expString,t,c-1,ru-1)),10);
            mycell.position = PCsomata(c,1:3);
%             mycell.clusterID = Clusters_str(c,connProbSelect);
            PCcells_str{c,ru}=mycell;
%             PCcells_str_spk{c,ru} = load(sprintf('%sexperiment_%d/%s/%d_%d_%d_st.txt',mypath,aa,expString,t,c-1,ru-1));
%             PCcells_str_dist{c,ru} = load(sprintf('%sexperiment_%d/%s/%d_%d_%d_DIST.txt',mypath,aa,expString,t,c-1,ru-1));
        end
        
        for c=pc+1:pv
            mycell = ncell(load(sprintf('%sexperiment_%d/%s/%d_%d_%d.txt',mypath,aa,expString,t,c-1,ru-1)),10);
            PVcells_str{c,ru}=mycell;%.hasPersistent(1000,8,mycell.tstop-1001); % paper?
        end
        for c=pv+1:cb
            mycell = ncell(load(sprintf('%sexperiment_%d/%s/%d_%d_%d.txt',mypath,aa,expString,t,c-1,ru-1)),10);
            CBcells_str{c,ru}=mycell;%.hasPersistent(1000,8,mycell.tstop-1001); % paper?
        end
        for c=cb+1:cr
            mycell = ncell(load(sprintf('%sexperiment_%d/%s/%d_%d_%d.txt',mypath,aa,expString,t,c-1,ru-1)),10);
           CRcells_str{c,ru}=mycell;%.hasPersistent(1000,8,mycell.tstop-1001); % paper?
        end
    end
    RUNS_str{1,t} = PCcells_str(:,:);
end

runtime = toc
fprintf('DONE!\n');

%% Plot
% close all;

ru = 1 ;% what run of this cluster

FFS = [];
FFPA = [];
for t=1:max(Clusters_str(1:nPC,stateID))
    for c=1:nPC
%         comprobj(c) = RUNS_str{1,t}{c,ru};
        FFS(t,c) = RUNS_str{1,t}{c,ru}.spike_count(500,1500) ;
        FFPA(t,c) = RUNS_str{1,t}{c,ru}.spike_count(1500,tstop) / 3.5;
    end
end
% figure();
% plot(FFS');hold on;
figure();
plot(FFPA');hold on;
plot(std(FFPA),'k','linewidth',2);hold on;

% tmp=find(StimVect_str);
% for c = 1:length(tmp)
% plot(tmp(c),RUNS_str{1,t}{tmp(c),ru}.spike_count(1500,tstop) / 3.5,'+r');hold on;
% end
%% Firing Frequency:
close all;
PCFF_mean = [];

for dataset=1
    % load dataset:
    if dataset == 1
        NC = max(Clusters_str(1:nPC,stateID));
        RUNS = RUNS_str;
    else 
        NC = NC_rnd(Sid);
%         NC = NC_str(Sid);
        RUNS = RUNS_rnd;
    end



%     % if in random network, replace the clusterID with that of the
%     % structured, becuse we dont have real clusters in random condition
%     if dataset == 2
%        for clu = 1:NC_str
%            for ru = 1:exprun % what run of this cluster
%                 for c=1:pc
%                     RUNS{1,clu}{c,ru}.clusterID = RUNS_str{1,clu}{c,ru}.clusterID;
%                 end
%            end
%        end
%     end
    
    
    
    FF = [];
    TMP=[];
    for ru = exprun % what run of this cluster
        FFPerClust=cell(NC(Sid));
        PAPerClust=zeros(NC(Sid));

        for stimclu =1:NC(Sid)
            for cl=1:NC(Sid)
                for c=1:pc
                    if(RUNS{1,stimclu}{c,ru}.clusterID == cl)
%                         if(RUNS{1,stimclu}{c,ru}.persistentActivity)
                            %                         PAPerClust(cl,stimclu) = PAPerClust(cl,stimclu) +1;
                            responce = RUNS{1,stimclu}{c,ru}.mv((stimend*10)+1:end);
                            FFPerClust{cl,stimclu} = [FFPerClust{cl,stimclu} ,length(find( ([0;diff(sign(diff(responce)))<0;0] & [sign(responce)==1]) )) / ((RUNS{1,stimclu}{c,ru}.tstop-(stimend+1))/1000) ];
%                         end
                    end
                end
            end
        end
        %     PAP(:,:,ru) = (PAPerClust.*(~eye(length(PAPerClust))));
        %     FFPC(:,:,ru) = cellfun(@(x,i) i*(nanmean(x)),FFPerClust,mat2cell(cellfun(@(x) ~isempty(x), FFPerClust),ones(NC(Sid),1),ones(NC(Sid),1)) );


        TMP = cellfun(@(x) nanmean(x),FFPerClust );
        TMP(~isfinite(TMP))  = 0;
%         FF(:,:,ru) = TMP./ ((RUNS{1,stimclu}{c,ru}.tstop-(stimend+1))/1000)  ; % Get Freq, no num of spikes
        FF(:,:,ru) = TMP;
    end

    FFnonstim = FF;

    FFstim = [];
    for i=exprun
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

% figure('Name','Random Firing Frequency' );
% imagesc(PCFF_mean(:,:,2));
% colormap(ccm(1:round(max(max(PCFF_mean(:,:,2)))*100),:));
% colorbar;
% set(gcf,'Color',[1,1,1]);
% 
% figure('Name','Histo Firing Frequency' );
% ph = bar( 0:0.5:ceil(max(max(max(PCFF_mean)))), histc(reshape(PCFF_mean(:,:,1),1,[]),0:0.5:ceil(max(max(max(PCFF_mean)))) )' );hold on;
% ch = get(ph,'Children');
% set(ch,'FaceAlpha',0.6, 'FaceColor', 'b') ;
% ph = bar( 0:0.5:ceil(max(max(max(PCFF_mean)))), histc(reshape(PCFF_mean(:,:,2),1,[]),0:0.5:ceil(max(max(max(PCFF_mean)))) )' );hold on;
% ch = get(ph,'Children');
% set(ch,'FaceAlpha',0.6, 'FaceColor', 'r') ;
% set(gcf,'Color',[1,1,1]);
% 
% figure;
% bar(cpc);
% set(gcf,'Color',[1,1,1]);
%% load morphologies:
% addpath('~/Documents/MATLAB/Dendrites/');
% addpath(genpath('~/Documents/MATLAB/Dendrites/TREES1.15'));
% morphPath = '/srv/userdata/HOMES/stefanos/Documents/MATLAB/Dendrites/' ;
% 
% % Load morphologies names from Smith's lab
% morphDir = 'SmithMorphologies';
% names = dir([morphPath,morphDir,'/']);
% names = {names(~[names.isdir]).name};
% names = names(cellfun(@(x) strcmp(x(end-3:end),'.swc'), names));
% 
% % Load all morphologies as ncell objects:
% for i=1:length(names)
%     i
%     PC(i) =  ncell([],1,[morphPath,morphDir,'/',names{i}]);
% end


%% Plot
% close all;

ru = 1 ;% what run of this cluster
Incommers = {};
for c=1:nPC
    Incommers{c} = find(AllConnMat_str(1:nPC,c));
end
tmpArr=[];
tmpMax=[];
tmpEachArr={};
for c=1:nPC
    tmpEachArr{c}=[];
    for i=1:length(Incommers{c})
        tmpArr = [tmpArr,PCcells_str_dist{Incommers{c}(i),ru}];
        tmpEachArr{c} = [tmpEachArr{c},PCcells_str_dist{Incommers{c}(i),ru}];
    end
    tmpMax(c) = max(tmpEachArr{c});
end
% tmpMin = min(tmpArr);
% tmpMax = max(tmpArr);
cm = jet( ceil(max(tmpArr)) - floor(min(tmpArr)) );
cm(1,:) = [1,1,1];
gm = gray(200);
gm = gm(end:-1:100,:); % overhead; reversed: not starting from black
colmap = [cm;gm];

figure();
for c=1:nPC
    colormap(colmap);
    c
% %     subplot('Position',[0.05,0.05,0.9,0.1]);
% subplot('Position',[0.05,0.55,0.9,0.4]);
% %     The raster plot below is VERY heavy...
% %     for syn = 1:length(Incommers{c})
% %         for l=1:size(PCcells_str_spk{Incommers{c}(syn),ru},1)
% %             plot([PCcells_str_spk{Incommers{c}(syn),ru}(l),PCcells_str_spk{Incommers{c}(syn),ru}(l)],[syn-1,syn],'k');hold on;
% %         end 
% %         HIST(syn) = PCcells_str_dist{Incommers{c}(syn),ru};
% %     end
%     imagesc( repmat( sum( incommingSynEvents{c} ) , length(Incommers{c}), [] ) + ( ceil(max(tmpArr)) - floor(min(tmpArr)) ) );
%     axis([0,tstop,0,length(Incommers{c})])
%     set(gca,'ytick',[]);
%     set(gca,'xtick',[]);
%     hold on;
    
    subplot('Position',[0.05,0.55,0.9,0.4]);
    plot(RUNS_str{1,t}{c,ru}.mv);
    axis([0,tstop*10,-90,60]);
    
    subplot('Position',[0.05,0.2,0.9,0.3]);
    incommingSynEvents{c}=zeros(length(Incommers{c}),tstop);
    
    for syn = 1:length(Incommers{c})
        for l=1:size(PCcells_str_spk{Incommers{c}(syn),ru},1)
            incommingSynEvents{c}(Incommers{c}(syn),round(PCcells_str_spk{Incommers{c}(syn),ru}(l))) = 1;
        end 
    end
    axis([0,tstop,0,length(Incommers{c})])
    ylabel('Incomming Synaptic Events');
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);
    imagesc(incommingSynEvents{c});
%     colormap( cm(1:round(tmpMax(c)) - floor(min(tmpArr)),:) );
    disp(sprintf('Mean distance from soma is: %f\n' , mean(tmpEachArr{c}) ))
    disp(sprintf('Incomming connections are: %f\n' , length(Incommers{c}) ))
    
    pause;
    cla;
end


%% detect UP states:

for c=1:nPC
    MUA(c,:) = RUNS_str{1,2}{c,1}.mv';
end
MUA = mean(MUA,1);

% get UP states more than 0.1 sec (for now):
[S,D,UPs]=findUPstates(MUA, 2, 3, RESTING, 1000 ) ;

% % Truncate all by the shorter UP state:
% UP_st = cell2mat( cellfun(@(x) x(1:min(cellfun(@length,UPs)))',UPs,'UniformOutput',false ) ) ;
% 
% UP_sa = mean(UP_st,1) ;
% 
% plot(UP_st', 'color', [0.9,0.9,0.9]);hold on;
% plot(UP_sa', 'k');


cl = lines(7);
IVu=[];
figure();
for t=2:8
    for i=1:length(D)
        VCi{i,t} = PCcells_str_i{t-1,1}.mv(S(i):S(i)+D(i));
    end
    % Truncate all by the shorter UP state:
    VCi_avr = cell2mat( cellfun(@(x) x(1:min(cellfun(@length,VCi(:,t))))',VCi(:,t),'UniformOutput',false ) ) ;
    VCi_avr = mean(VCi_avr,1) ;
    
    plot(VCi_avr, 'color',cl(t-1,:));hold on;
    plot([0,length(VCi_avr)],[min(VCi_avr),min(VCi_avr)], 'color',(cl(t-1,:)+0.5)/norm(cl(t-1,:)+0.5),'linewidth',2);hold on;
    lstr{t-1} = sprintf('Vclamp @%.1f',RESTING + test_vclamp(t-1));
    IVu(t-1,:) = min(VCi_avr) ./ RESTING + test_vclamp(t-1);
end
legend(lstr);




MUA_DOWN=[];
IVd=[];
figure();
for t=2:8
    MUA_DOWN(t-1,:) = PCcells_str_i{t-1,1}.mv(21350:24780);
    plot(MUA_DOWN(t-1,:), 'color',cl(t-1,:));hold on;
    plot([0,length(MUA_DOWN(t-1,:))],[min(MUA_DOWN(t-1,:)),min(MUA_DOWN(t-1,:))], 'color',(cl(t-1,:)+0.5)/norm(cl(t-1,:)+0.5),'linewidth',2);hold on;
    IVd(t-1,:) = min(MUA_DOWN(t-1,:)) / RESTING + test_vclamp(t-1);
end
legend(lstr);


figure();
plot(IVu);hold on;
max(diff(IVu))
plot(IVd,'r');hold on;
max(diff(IVd))
%%

plot( PCcells_str{c,ru}.mv & ((PCcells_str{c,ru}.mv-RESTING)>5) );
hold on;
plot(PCcells_str{c,ru}.mv, 'r')

%%

t=1; % Initial without vclamps
VCLAMP = 0;
GPYID = 0;

fprintf('Running NEURON... ');
unix(sprintf('kill -9 `ps -ef | grep stefanos | grep nrniv | grep -v grep | awk ''{print$2}''`'));
tic;
%         [nrnStatus,nrnCmdOut{t}] = unix('./parallel_pfc');
EXPERIMENT = 1;  % 0 = random, 1= structured gia to $#$#% NEURON
%         [nrnStatus,nrnCmdOut{t}] = unix(sprintf('mpiexec -np 8 %smechanism_complex/x86_64/special -mpi -nobanner -c "PARALLEL=1" -c "SIMPLIFIED=0" -c "CLUSTER_ID=%d" -c "EXPERIMENT=%d"  -c "AA=%d" final.hoc',mypath2,t,EXPERIMENT,aa));
[nrnStatus,nrnCmdOut{t}] = unix(sprintf('mpirun -n 24 -mca btl ^openib %smechanism_complex/x86_64/special -mpi -nobanner -c "PARALLEL=1" -c "SIMPLIFIED=0" -c "CLUSTER_ID=%d" -c "EXPERIMENT=%d"  -c "AA=%d" -c "VCLAMP=%f"  -c "GPYID=%d" final.hoc',mypath2,t,EXPERIMENT,aa, VCLAMP, GPYID));
runtime = toc



%% detect UP states; Current clamp:

for c=2:nPC
    MUA(c,:) = RUNS_str{1,1}{c,1}.mv';
end
MUA = mean(MUA,1);

% get UP states more than 0.1 sec (for now):
[S,D,UPs]=findUPstates(MUA, 2, 3, RESTING, 1000 ) ;


for i=1:length(D)
    VC{i} = RUNS_str{1,t}{1,1}.mv(S(i):S(i)+D(i));
end
% Truncate all by the shorter UP state:
VC_avr = cell2mat( cellfun(@(x) x(1:min(cellfun(@length,VC)))',VC,'UniformOutput',false )' ) ;
VC_avr = mean(VC_avr,1) ;

plot(VC_avr)

Rin_active = (VC_avr(1) - min(VC_avr)) / (-0.2) ;

