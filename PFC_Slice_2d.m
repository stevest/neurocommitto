%%  ---- Initialize Pyramidal cells ---- 
nPC = 65;%216;

% Initialize somata positions
% PCsomata = CreateSliceNetwork([4,4,4], 100, 30);
PCsomata = CreateRandomNetwork(nPC, 200, 3);
% scatter3(PCsomata(:,1),PCsomata(:,2), PCsomata(:,3), 20, [0,0,1], 'fill', 'o');
% axis equal; hold on;

distPC2PC = generateDistanceMat(PCsomata', 0);
connBinsPC2PC = 20:30:500;
connProbsPC2PC = 0.25 .* exp(-0.006 * connBinsPC2PC);

recipBinsPC2PC = 20:30:500;
recipProbsPC2PC = 0.12 .* exp(-0.006 * recipBinsPC2PC);
% Update probabilities from Perin et al. (Figure S4)
NconnBins = [0,1,2,3];
NincomingProbs = [0.1, 0.2, 0.25, 0.4] ;
NoutgoingProbs = [0.12, 0.25, 0.24, 0.2] ;

% % Yuste Fig4D:
% bar(histc(reshape(ConnMatPV2PC.*distPV2PC,1,[]),1:20:600));

% Pyramidals connect to all
glProbConn = 0.7;
PC2PC = zeros(nPC,nPC,4);
PC2PC(:,:,5) = create_graph_WS(distPC2PC,glProbConn,0.99); % 1.0=Random graph, 0.0=Watts-Strogatz graph
PC2PC(:,:,4) = create_graph_WS(distPC2PC,glProbConn,0.5); % 1.0=Random graph, 0.0=Watts-Strogatz graph
PC2PC(:,:,3) = create_graph_WS(distPC2PC,glProbConn,0.01); % 1.0=Random graph, 0.0=Watts-Strogatz graph
PC2PC(:,:,2) = create_graph_DD(distPC2PC,glProbConn,connBinsPC2PC, connProbsPC2PC);
PC2PC(:,:,1) = create_graph_CN(distPC2PC,glProbConn, connBinsPC2PC, connProbsPC2PC, ...
    recipBinsPC2PC, recipProbsPC2PC,NconnBins,NincomingProbs,NoutgoingProbs);

% clustering = 3;
for i=1:size(PC2PC,3)
    tempPC2PC(:,:,i) = PC2PC(:,:,i);
end

%% ---- affinity propagation similarity measure: ---- 
for d=1:size(tempPC2PC,3)
    % Find nearest neighbors of Pyramidals only
    [CNi,CNo] = commonNeighbors(tempPC2PC(:,:,d));
    onlyCN{d} = [CNi + CNo] .* tempPC2PC(:,:,d);
    
    %perform affinity propagation algorithm:
    [NC{d},labels{d},Sid{d}] = runAffinityPropagation(onlyCN{d});
    
    %how many cells in each cluster?
    cellsPerCluster(d) = {histc(labels{d}(:,Sid{d}),1:NC{d}(Sid{d}))'};
%     bar(cellsPerCluster{d})

    
    Deg{d} = [];
    for i=1:NC{d}(Sid{d})
        tempIdx = find(labels{d}(:,Sid{d}) == i);
        tempMat = tempPC2PC(tempIdx,tempIdx);
%         CC(i) = clust_coeff(tempMat);
        Deg{i,d} = [Deg{d},degrees(tempMat)];
    end
    [~,IDX(d)]=max(cellfun(@mean,Deg(:,d)));
    % [~,IDX]=max(CC);
    
    % Get cluster with greater degree:
    [r,~,~] = find(labels{d}(:,Sid{d}) == IDX(d));
    StimVect(d) = {zeros(nPC,1)};
    StimVect{d}(r) = 1;
end
%%  ---- Isolate Cluster ---- 
equalizeClusters(cellsPerCluster,thrshld)

for d=1:size(tempPC2PC,3)
    exportPC2PC(d) = {tempPC2PC(find(StimVect{d}),find(StimVect{d}),d)};
    exportPCsomata(d) = {PCsomata(find(StimVect{d}),:)};
    exportStimVect(d) = {ones(size(r,1),1)};
    nPC = length(PCsomata);
end
%%  ---- Rest of the cells ---- 
% tempPC2PC = PC2PC(:,:,4); % replace connectivity to your liking
% tempPC2PC(r,r) = PC2PC(r,r,1);
clstr = 1;
% Rest of the cells!
nPV = round(nPC*17/65);
nCB = round(nPC*9/65);
nCR = round(nPC*9/65);
nAllCells = sum([nPC,nPV,nCB,nCR]);
AllCells = [nPC , nPV , nCB , nCR, nAllCells];


PVsomata = CreateCubeNetworkPV(250, nPV); % 226.6 per mm squared (?!?) paper?
% scatter3(PVsomata(:,1),PVsomata(:,2), PVsomata(:,3), 30, [0.5,0.5,0.5], 'fill', 'o');
% axis equal;

CBsomata = CreateCubeNetworkPV(200, nCB);
CRsomata = CreateCubeNetworkPV(200, nCR);

% Distance of each interneuron from Pyramidal and connect based on
% probability from Yuste 2011:
distPV2PC = distancePV2PC(PVsomata,PCsomata);
connBinsPV2PC = 0:20:500;
connProbsPV2PC = linspace(1,0,length(connBinsPV2PC));
ConnMatPV2PC = connectPV2PC(distPV2PC,connBinsPV2PC,connProbsPV2PC);

% Populate the final all-to-all connectivity matrix:
AllConnMat = zeros(nAllCells);
% Set indies for ease of mind:
pc = size(PCsomata,1);
pv = size(PCsomata,1) + size(PVsomata,1);
cb = size(PCsomata,1) + size(PVsomata,1) + size(CBsomata,1);
cr = size(PCsomata,1) + size(PVsomata,1) + size(CBsomata,1) + size(CRsomata,1);

% tempPC2PC(sub2ind([nPC,nPC],[1:pc],[1:pc])) = 1; % autapses in Pyramidals
tmptmp = tempPC2PC(:,:,clstr);
tmptmp(sub2ind([nPC,nPC],[1:pc],[1:pc])) = 1;
AllConnMat(1:pc,1:pc) = tmptmp;


% Pyramidals to all types of interneurons:
AllConnMat(1:pc,pc+1:end) = 1; % Connected to all
% PVs connect to all other PVs + autapses (need to gap junctions?)
AllConnMat(pc+1:pv,pc+1:pv) = 1;
% PVs to PCs based on above connectivity (Yuste 2011):
AllConnMat(pc+1:pv,1:pc) = ConnMatPV2PC;
% CB only connect to PC (Xenia) na to psa3w...
AllConnMat(pv+1:cb,1:pc) = 1;
% CR only connect to PC and CB (Xenia) na to psa3w...
AllConnMat(cb+1:cr,1:pc) = 1;
AllConnMat(cb+1:cr,pv+1:cb) = 1;

figure();
imagesc(AllConnMat)
% clust_coeff(PC2PC(:,:,4))

save('../CRUN.mat');

%%  ---- Export parameters to NEURON ---- 
mypath = '/srv/userdata/HOMES/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/';
exportNetworkParameters(AllCells,AllConnMat,mypath);
exportStimulationParameters(AllCells,StimVect{clstr},mypath);
%%  ---- Plot clustering ---- 
figure('renderer','opengl');hold on;
for i=1:NC{clstr}(Sid{clstr})
    clustCol = rand(1,3);
    tempIdx = find(labels{clstr}(:,Sid{clstr}) == i);
    for c1 = 1:length(tempIdx)
        for c2 = 1:length(tempIdx)
            if AllConnMat(tempIdx(c1),tempIdx(c2)) && (i == IDX(clstr))
                plot3( PCsomata([tempIdx(c1),tempIdx(c2)],1),...
                    PCsomata([tempIdx(c1),tempIdx(c2)],2),...
                    PCsomata([tempIdx(c1),tempIdx(c2)],3), 'color', clustCol );
            end
        end
    end
    scatter3(PCsomata(tempIdx,1),PCsomata(tempIdx,2), PCsomata(tempIdx,3),...
        45, clustCol, 'fill', 'o');hold on;
end
axis equal; hold on;
%%  ---- Perin's statistics ---- 
% Use custom combinations function for effficiency:
% idxCombs3 = [combntns(1:3,2); combntns(3:-1:1,2)];
% cluster3= [];
% dist3 = [];
idxCombs6 = [combntns(1:6,2); combntns(6:-1:1,2)];
cluster6= [];
dist6 = [];

%Cluster of 6 pyramidals
ctr=1;
for t = 1:30
    layerDepth = rand(1) * 250;
    dozen = find((PCsomata(:,3)>(layerDepth-1)) & (PCsomata(:,3)<(layerDepth+1)));
    
    if(length(dozen) >= 6)
        allCombs6 = combntns(1:length(dozen),6);
        for i=1:size(allCombs6,1)
            counter = 0;
            meanDist = [];
            for j=1:length(idxCombs6)
                tempIdx = allCombs6(i,idxCombs6(j,:));
                
                if AllConnMat( dozen(tempIdx(1)), dozen(tempIdx(2)) )
                    counter = counter + 1;
                    meanDist = [meanDist,distPC2PC( dozen(tempIdx(1)), dozen(tempIdx(2)) )] ;
                end
            end
            cluster6(ctr) = counter;
            dist6(ctr) = mean(meanDist);
            ctr = ctr+1;
        end
    end
    
%     if(length(dozen) >= 3)
%         allCombs3 = combntns(1:length(dozen),3);
%         for i=1:size(allCombs3,1)
%             counter = 0;
%             meanDist = [];
%             for j=1:length(idxCombs3)
%                 tempIdx = allCombs3(i,idxCombs3(j,:));
%                 
%                 if AllConnMat( dozen(tempIdx(1)), dozen(tempIdx(2)) )
%                     counter = counter + 1;
% %                     meanDist = [meanDist,distPC2PC( dozen(tempIdx(1)), dozen(tempIdx(2)) )] ;
%                 end
%             end
%             cluster3(ctr) = counter;
% %             dist3(ctr) = mean(meanDist);
%             ctr = ctr+1;
%         end
%     end
end

% Initialize plot:
PerinPlot={};
PerinPlot(1,1:length(50:25:250))={[]};


for i=1:ctr-1
    b = find( histc(dist6(i), 50:25:250) );
    if ~isempty(b)
    PerinPlot{1,b} = [PerinPlot{1,b}, cluster6(i)] ;
    end
end

bar(cellfun(@nanmean,PerinPlot));


% [aa,b] = hist(distPC2PC(:), length(idxCombs6)+1);
% aa=aa/max(aa);
% figure;plot(mean(clusterHist6)./aa);figure(gcf);
%% Temp tests
clear DEG_R DEG_C CC_R CC_C CC_D DEG_D CC_CN DEG_CN
tempConnMat = [];
for i=1:100
    i
    %     Watts & Strogazt random
    tempConnMat = create_graph_WS(distPC2PC,0.99); % 1.0=Random graph, 0.0=Watts-Strogatz graph
    CC_R(i,:) = clust_coeff(tempConnMat);
    DEG_R(i,:) = sort(degrees(tempConnMat),'descend');
    %     Watts & Strogazt clustered
    tempConnMat = create_graph_WS(distPC2PC,0.01); % 1.0=Random graph, 0.0=Watts-Strogatz graph
    CC_C(i,:) = clust_coeff(tempConnMat);
    DEG_C(i,:) = sort(degrees(tempConnMat),'descend');
    %     Distance dependent:
    tempConnMat = create_graph_DD(distPC2PC, connBinsPC2PC, connProbsPC2PC);
    CC_D(i,:) = clust_coeff(tempConnMat);
    DEG_D(i,:) = sort(degrees(tempConnMat),'descend');
    %     Common neighbour dependent (Perin et al)
    tempConnMat = create_graph_CN(distPC2PC, connBinsPC2PC, connProbsPC2PC, ...
        recipBinsPC2PC, recipProbsPC2PC,NconnBins,NincomingProbs,NoutgoingProbs);
    CC_CN(i,:) = clust_coeff(tempConnMat);
    DEG_CN(i,:) = sort(degrees(tempConnMat),'descend');
end

plot(mean(DEG_R),'g');hold on;
plot(mean(DEG_C),'r');hold on;
plot(mean(DEG_D),'m');hold on;
plot(mean(DEG_CN),'k');hold on;

mean(CC_R)
mean(CC_C)
mean(CC_D)
mean(CC_CN)



%%
DistMat = generateDistanceMat(PCsomata', 0);

% Import initial connectivity probabilities from Perin et al. (figure 1)
connBins = 20:30:400;
connProbs = 0.25 .* exp(-0.006 * connBins);
recipBins = 20:30:400;
recipProbs = 0.12 .* exp(-0.006 * recipBins);
% connBins = [20,50,90,125,160,190,225,260,295,340];
% connProbs = [0.25, 0.21, 0.18, 0.15, 0.11, 0.09, 0.06, 0.04, 0.02, 0.01];
% recipBins = [20,50,90,125,160,190,225,260,295,340];
% recipProbs = [0.12, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.003];


% connProbs = connProbs/sum(connProbs);
% recipProbs = recipProbs/sum(recipProbs);

% Update probabilities from Perin et al. (Figure S4)
NconnBins = [0,1,2,3];
NincomingProbs = [0.1, 0.2, 0.25, 0.4] ;
NoutgoingProbs = [0.12, 0.25, 0.24, 0.2] ;

% Initialize network: Create connectivity based on intersomatic distance:
% Histo of distances:
% [a,b] = hist(DistMat(:), connBins);
% a=a/max(a);
a=ones(1,length(connBins));
ConnMat = initializeNetwork(DistMat',a,connBins,connProbs,recipBins,recipProbs);
sc(ConnMat);

tic;
ConnMat = cu_reshapeNetwork(DistMat',ConnMat',0.1,NconnBins,NincomingProbs, NoutgoingProbs,recipBins,recipProbs, 1); % Probabilities map(probsMat) is optional
toc

% tic;
% [ConnMat, probsMat, commNeighs] = reshapeNetwork(DistMat',ConnMat',a,NconnBins,NincomingProbs, NoutgoingProbs,recipBins,recipProbs, 0.2); % Probabilities map(probsMat) is optional
% toc
% sc(ConnMat)

figS7 = histc(commNeighs(find(commNeighs)),0:1:300);
figS7 = figS7/sum(figS7);
bar(figS7)
[~,~,CC] = clust_coeff(ConnMat) ;

ConnMat = zeros(size(Points, 1));
ConnMat(round(rand(120835,1)*4826808)) = 1;

mean(sum(ConnMat))

% check connectivity fadeout:
for i=1:3446 % Perin
    SimpleConns(i) = ConnMat(ceil(rand(1)*length(Points)),ceil(rand(1)*length(Points))) ;
    SimpleDists(i) = DistMat(ceil(rand(1)*length(Points)),ceil(rand(1)*length(Points))) ;
end
SingleConn = histc( SimpleDists(find(SimpleConns)), connBins );
SingleConn = SingleConn / sum(SingleConn);
bar(SingleConn);


% Use custom combinations function for effficiency:
sixOfTwelve = combntns(1:12,6);
idxCombs6 = [combntns(1:6,2); combntns(6:-1:1,2)];

% Mou kanei bias: sample apo liga mono kyttara!!!! (ane3artita apo to
% posous syndiasmous exoun..)
dozens = ceil(rand(100,12) * length(Points));
% dozens = combntns(randomSample,12);

%Motifs of 6 pyramidals
mot6=[];
motBins = [50, 75, 100, 125, 150, 175, 200, 225];
for m = 1:length(motBins)
    motif6 = zeros(length(dozens),length(sixOfTwelve) );
    for t = 1:length(dozens)
        for i=1:length(sixOfTwelve)
            avgDist = 0;
            for j=1:length(idxCombs6)
                tempIdx = sixOfTwelve(i,idxCombs6(j,:));
                avgDist = avgDist + DistMat( dozens(t,tempIdx(1)), dozens(t,tempIdx(2)) );
            end
            avgDist = avgDist / length(idxCombs6);
            if((avgDist>(motBins(m)-25)) && (avgDist<motBins(m)))
                counter=0;
                for j=1:length(idxCombs6)
                    tempIdx = sixOfTwelve(i,idxCombs6(j,:));
                    counter = counter + ConnMat( dozens(t,tempIdx(1)), dozens(t,tempIdx(2)) );
                end
                motif6(t,i) = counter;
            end
        end
    end
    mot6(m) = mean(motif6(find(motif6>0)));
end
bar(mot6)


%observed connections:
twelveOfTwelve = combntns(1:12,12);
idxCombs12 = [combntns(1:12,2); combntns(12:-1:1,2)];
mot6=[];

motif6 = zeros(length(dozens),length(sixOfTwelve) );
for t = 1:length(dozens)
    for i=1:length(sixOfTwelve)
        counter=0;
        for j=1:length(idxCombs6)
            tempIdx = sixOfTwelve(i,idxCombs6(j,:));
            counter = counter + ConnMat( dozens(t,tempIdx(1)), dozens(t,tempIdx(2)) );
        end
        motif6(t,i) = counter;
    end
end
mot6 = histc(motif6,0:1:12);
mean(mot6')


%% Clustering in larger groups
% Use custom combinations function for effficiency:
allCombs3 = combntns(1:12,3);
allCombs6 = combntns(1:12,6);
allCombs8 = combntns(1:12,8);
idxCombs3 = [combntns(1:3,2); combntns(3:-1:1,2)];
idxCombs6 = [combntns(1:6,2); combntns(6:-1:1,2)];
idxCombs8 = [combntns(1:8,2); combntns(8:-1:1,2)];


%Cluster of 6 pyramidals
for t = 1:length(dozens)
    for i=1:length(allCombs6)
        counter = 0;
        for j=1:length(idxCombs6)
            tempIdx = allCombs6(i,idxCombs6(j,:));
            counter = counter + ConnMat( dozens(t,tempIdx(1)), dozens(t,tempIdx(2)) );
        end
        cluster6(t,i) = counter;
    end
    clusterHist6(t,:) = histc(cluster6(t,:),0:length(idxCombs6)) ;
    clusterHist6(t,:) = clusterHist6(t,:)/norm(clusterHist6(t,:));
end
[aa,b] = hist(DistMat(:), length(idxCombs6)+1);
aa=aa/max(aa);
figure;plot(mean(clusterHist6)./aa);figure(gcf);
%% Call Affinity Propagation (on common neighbors as similarity):
PCsomata = [CreateRandomNetwork(40, 10, 3);CreateRandomNetwork(40, 40, 3);CreateRandomNetwork(40, 80, 3)];

% programs for adaptive Affinity Propagation clustering; an improvement
% of Affinity Propagation clusteirng (see Frey & Dueck, Science, Feb. 2007)
% Note: Statistics Toolbox of Matlab needs to be installed
% WANG Kaijun: wangkjun@yahoo.com, July, Sept. 2007.

algorithm = 1;  % 1 --- adaptive AP, 0 --- original AP
nrun = 5000;   % max iteration times, default 50000
nrun2 = 200;   % max iteration times for original AP
nconv = 50;     % convergence condition, default 50
pstep = 0.01;   % decreasing step of preferences: pstep*pmedian, default 0.01
lam = 0.5;      % damping factor, default 0.5
cut = 3;        % after clustering, drop an cluster with number of samples < cut
% splot = 'plot'; % observing a clustering process when it is on
splot = 'noplot';
M=[];
p = [];
type = 1;
simatrix=0;
data = PCsomata;
[nrow, dim] = size(data);
% taking true class labels from a data file
truelabels = ones(nrow,1);


tic;

[labels,NCs,labelid,iend,Sp,Slam,NCfixs] = adapt_apcluster(data,type,...
    p,pstep,simatrix,'convits',nconv,'maxits',nrun,'dampfact',lam,splot);

[NC,Sil,Silmin] = solution_evaluation(data,M,labels,NCs,...
    NCfixs,simatrix,nrow,type,cut);
trun = toc;

fprintf('\n## Running time = %g seconds \n', trun);
fprintf('## Running iterations = %g \n', iend);

% finding an optimal clustering solution
[Smax, Sid] = max(Sil);
NCopt = NC(Sid);
fprintf('\n## Clustering solution by adaptive Affinity Propagation:\n');
fprintf('  Optimal number of clusters is %d, Silhouette = %g,\n',NCopt,Smax);
if Smax < 0.3
    R = length(NC);
    R = ceil(R/2):R;
    [Tmax, Q] = max(Silmin(R));
    Sid = R(Q);
    NCopt2 = NC(Sid);
    fprintf('  If Silhouette values are small & NCs are large, Optimal NC is %d,\n',NCopt2);
    fprintf('  where min Silhouette of single cluster is %g.\n',Tmax);
    fprintf('  The optimal solution (class labels) is in labels(:,Sid)');
end
fprintf('\n## Silhouette values at different NCs: [NC;Sil;Silmin] \n');
disp([NC;Sil;Silmin]);



truek = unique(truelabels);
truek = length(truek);
if truek > 1
    C = valid_external(labels(:,Sid), truelabels);
    fprintf('Fowlkes-Mallows validity index: %f\n', C(4));
end
if NCopt == truek
    fprintf('\n## Error rate of clustering solution might be inaccurate if large');
    fprintf('\n     (then use FM index instead) and it is for reference only:');
    valid_errorate(labels(:,Sid), truelabels);
end
%clf; plotdata_bylabels(data,truelabels,2,0,'co');
%clf; plotdata_bylabels(data,labels(:,M),2,0,'nb');

scatter3(PCsomata(:,1),PCsomata(:,2), PCsomata(:,3), 50, labels(:,find(NCs==NCopt)), 'fill', 'o');
axis equal; hold on;
