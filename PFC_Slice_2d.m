% Initialize somata positions 
% PCsomata = CreateSliceNetwork([4,4,4], 100, 30);
PCsomata = CreateRandomNetwork(65, 100, 3);
% scatter3(PCsomata(:,1),PCsomata(:,2), PCsomata(:,3), 5, [0,0,1], 'fill', 'o');
% axis equal; hold on;
PVsomata = CreateCubeNetworkPV(60, 17); % 226.6 per mm squared (?!?) paper?
% scatter3(PVsomata(:,1),PVsomata(:,2), PVsomata(:,3), 30, [0.5,0.5,0.5], 'fill', 'o');
% axis equal;

CBsomata = CreateCubeNetworkPV(200, 9);
CRsomata = CreateCubeNetworkPV(200, 9);

nPC = length(PCsomata);
nPV = length(PVsomata);
nCB = length(CBsomata);
nCR = length(CRsomata);
AllCells = [nPC , nPV , nCB , nCR];
nAllCells = sum(AllCells);

% Distance of each interneuron from Pyramidal and connect based on
% probability from Yuste 2011:
distPV2PC = distancePV2PC(PVsomata,PCsomata);
connBinsPV2PC = 0:20:500;
connProbsPV2PC = linspace(1,0,length(connBinsPV2PC));
ConnMatPV2PC = connectPV2PC(distPV2PC,connBinsPV2PC,connProbsPV2PC);

distPC2PC = generateDistanceMat(PCsomata', 0);
% % Yuste Fig4D:
% bar(histc(reshape(ConnMatPV2PC.*distPV2PC,1,[]),1:20:600));

% Populate the final all-to-all connectivity matrix:
AllConnMat = zeros(nAllCells);
% Set indies for ease of mind:
pc = length(PCsomata);
pv = length(PCsomata) + length(PVsomata);
cb = length(PCsomata) + length(PVsomata) + length(CBsomata);
cr = length(PCsomata) + length(PVsomata) + length(CBsomata) + length(CRsomata);
%Pyramidals connect to all
AllConnMat(1:pc,1:pc) = create_graph(distPC2PC,0.9); % 0=Random graph, 4=Watts-Strogatz graph 
%Pyramidals to all types of interneurons:
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

sc(AllConnMat)

%% Export parameter matrices in .hoc file:

% Write NMDA results to a .hoc file:
fid = fopen('../importNetworkParameters.hoc','w');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');
fprintf(fid,'// Override variables\n');

fprintf(fid,sprintf('nPCcells=%d\n',nPC));
fprintf(fid,sprintf('nPVcells=%d\n',nPV));
fprintf(fid,sprintf('nCBcells=%d\n',nCB));
fprintf(fid,sprintf('nCRcells=%d\n',nCR));
fprintf(fid,sprintf('nAllCells=%d\n\n',nAllCells));

fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref connMatrix, weightsMatrix\n');
fprintf(fid,'connMatrix = new Matrix(nAllCells, nAllCells)\n');
fprintf(fid,'weightsMatrix = new Matrix(nAllCells, nAllCells)\n');

fprintf(fid,'\n\n// Import parameters: (long-long text following!)\n\n');
% network connectivity:
for i=1:length(AllConnMat)
    for j=1:length(AllConnMat)
        fprintf(fid,'connMatrix.x[%d][%d] = %d\n',i-1,j-1, AllConnMat(i,j));
    end
end
% Network synaptic weights
for i=1:length(AllConnMat)
    for j=1:length(AllConnMat)
        fprintf(fid,'weightsMatrix.x[%d][%d] = %f\n',i-1,j-1, AllConnMat(i,j));
    end
end
fclose(fid);
%%
% affinity propagation similarity measure:

% Find nearest neighbors of Pyramidals only
NN = nearestNeighbors(AllConnMat(1:pc,1:pc));
onlyNN = NN.*AllConnMat(1:pc,1:pc);
% NN(find(NN.*AllConnMat(1:pc,1:pc)));

%perform affinity propagation algorithm:

% Affinity Propagation clusteirng (see Frey & Dueck, Science, Feb. 2007)
% Note: Statistics Toolbox of Matlab needs to be installed

algorithm = 1;  % 1 --- adaptive AP, 0 --- original AP
nrun = 50000;   % max iteration times, default 50000
% nrun2 = 2000;   % max iteration times for original AP
nconv = 50;     % convergence condition, default 50
pstep = 0.01;   % decreasing step of preferences: pstep*pmedian, default 0.01
lam = 0.5;      % damping factor, default 0.5
cut = 3;        % after clustering, drop an cluster with number of samples < cut
%splot = 'plot'; % observing a clustering process when it is on
splot = 'noplot';

% initialization
type = 1;       % 1: Euclidean distances
simatrix = 1;   % 0: data as input; 1: similarity matrix as input
p = [];
Ms = [];
% derive similarity matrix
[rows,cols]=find(onlyNN);
for i=1:length(rows)
    M(i,1) = rows(i);
    M(i,2) = cols(i);
    M(i,3) = onlyNN(rows(i),cols(i));
end
data = [];
truelabels = ones(pc,1);


disp(' '); disp(['==> Clustering is running, please wait ...']);
tic;
[labels,NCs,labelid,iend,Sp,Slam,NCfixs] = adapt_apcluster(M,type,...
p,pstep,simatrix,'convits',nconv,'maxits',nrun,'dampfact',lam,splot);

[NC,Sil,Silmin] = solution_evaluation(data,M,labels,NCs,...
NCfixs,simatrix,pc,type,cut);
trun = toc;

fprintf('\n## Running time = %g seconds \n', trun);
fprintf('## Running iterations = %g \n', iend);

% finding an optimal clustering solution
solution_findK

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

for i=1:NC(Sid)
    tempIdx = find(labels(:,Sid) == i);
    tempMat = AllConnMat(tempIdx,1:pc);
    tempMat = tempMat(:,tempIdx);
    clust_coeff(tempMat)
    scatter3(PCsomata(tempIdx,1),PCsomata(tempIdx,2), PCsomata(tempIdx,3), 45, rand(1,3), 'fill', 'o');hold on;
end
axis equal; hold on;

clust_coeff(AllConnMat(1:pc,1:pc))
%% Export stimulation parameters in .hoc file:

% Write clustering results to a .hoc file:
fid = fopen('../importStimulationParameters.hoc','w');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');
fprintf(fid,'// Override variables\n');

fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref PcellStimList\n');
fprintf(fid,'PcellStimList = new Vector(nPCcells)\n');

fprintf(fid,'\n\n// Import parameters: (long-long text following!)\n\n');
% network stimulation:
for i=1:nPC
    for j=1:nPC
        fprintf(fid,'PcellStimList.x[%d] = %d\n',i-1,j-1, AllConnMat(i,j));
    end
end
fclose(fid);
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
