close all;clear all;
% Main parameters:
nPC = 14;  %was 14
nPV = 3; % was 6
nCB = 0;
nCR = 0;
AllCells = [nPC , nPV , nCB , nCR];
nAllCells = sum(AllCells);

% Create Network:
Points = CreateSobolNetwork(nAllCells, 28, 3); % was 28
% Points = CreateRandomNetwork(nAllCells, 320, 3); % was 28

figure('Renderer', 'OpenGL');
scatter3(Points(1:nPC,1), Points(1:nPC,2), Points(1:nPC,3),...
    '^','markerEdgeColor','b','markerFaceColor','b',...
    'SizeData',90); hold on;
scatter3(Points(nPC+1:nPC+nPV,1), Points(nPC+1:nPC+nPV,2), Points(nPC+1:nPC+nPV,3),...
    '^','markerEdgeColor','r','markerFaceColor','r',...
    'SizeData',50); hold on;
scatter3(Points(nPC+nPV+1:nPC+nPV+nCB,1), Points(nPC+nPV+1:nPC+nPV+nCB,2), Points(nPC+nPV+1:nPC+nPV+nCB,3),...
    '^','markerEdgeColor','r','markerFaceColor','r',...
    'SizeData',50); hold on;
scatter3(Points(nPC+nPV+nCB+1:nPC+nPV+nCB+nCR,1), Points(nPC+nPV+nCB+1:nPC+nPV+nCB+nCR,2), Points(nPC+nPV+nCB+1:nPC+nPV+nCB+nCR,3),...
    '^','markerEdgeColor','r','markerFaceColor','r',...
    'SizeData',50); hold on;
axis equal;

% Points = CreateRandomNetwork(pyramidalNo, 300, 3);
DistMat = generateDistanceMat(Points(1:nPC,:)');
mean(mean(DistMat))
max(max(DistMat))

% Use custom combinations function for effficiency:
% maxCombs = 10000;
allCombs3 = combntnsRND(nPC,3,4199);
allCombs4 = combntnsRND(nPC,4,8202);
allCombs5 = combntnsRND(nPC,5,11544);
allCombs6 = combntnsRND(nPC,6,12012);
allCombs7 = combntnsRND(nPC,7,9306);
allCombs8 = combntnsRND(nPC,8,5319);
% % 
idxCombs3 = [combntns(1:3,2); combntns(3:-1:1,2)];
idxCombs4 = [combntns(1:4,2); combntns(4:-1:1,2)];
idxCombs5 = [combntns(1:5,2); combntns(5:-1:1,2)];
idxCombs6 = [combntns(1:6,2); combntns(6:-1:1,2)];
idxCombs7 = [combntns(1:7,2); combntns(7:-1:1,2)];
idxCombs8 = [combntns(1:8,2); combntns(8:-1:1,2)];

% load('allCombs_30_6.mat');
% allCombs = combntns(1:pyramidalNo,3);
% allCombs3 = combntns(1:pyramidalNo,3);
% allCombs4 = combntns(1:nPC,4);
% allCombs5 = combntns(1:pyramidalNo,5);
% allCombs6 = combntns(1:pyramidalNo,6);
% allCombs7 = combntns(1:pyramidalNo,7);
% allCombs8 = combntns(1:pyramidalNo,8);
% idxCombs3 = [combntns(1:3,2); combntns(3:-1:1,2)];
% idxCombs4 = [combntns(1:4,2); combntns(4:-1:1,2)];
% idxCombs5 = [combntns(1:5,2); combntns(5:-1:1,2)];
% idxCombs6 = [combntns(1:6,2); combntns(6:-1:1,2)];
% idxCombs7 = [combntns(1:7,2); combntns(7:-1:1,2)];
% idxCombs8 = [combntns(1:8,2); combntns(8:-1:1,2)];

% From perin et al.
connBins = [20,50,90,125,160,190,225,260,295,340];
connProbs = [0.25, 0.21, 0.18, 0.15, 0.11, 0.09, 0.06, 0.04, 0.02, 0.01];
recipBins = [20,50,90,125,160,190,225,260,295,340];
recipProbs = [0.12, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.003];
% recipBins = [1,1000];
% recipProbs = [0.8, 0.8];

% NconnBins = [0,1,2,3,4];
% NconnProbs = [0.1, 0.21, 0.22, 0.28, 0.3];

NconnBins = [0,1,2,3];
NincomingProbs = [0.12, 0.21, 0.25, 0.4];
NoutgoingProbs = [0.1, 0.25, 0.22, 0.18];

%increased connectivity in Prefrontal
connProbs = connProbs * 1.3; % 10%-20% increase in connectivity
recipProbs = recipProbs * 1.8; % 80% increased bidirectional connectivity

% IR=0 ; % reciprocal probability is independend of connection probability!
ConnMat = initializeNetwork(DistMat',connBins,connProbs,recipBins,recipProbs);
disp(sprintf('Reciprocal connections: %3.2f%%\n',(sum(sum(ConnMat .* ConnMat')) * 100) / (nPC^2 - nPC) ));
disp(sprintf('Initial clustering coefficient (S) is: %f \n',nanmean(clustCoeff(ConnMat))));
[C1,C2,C]=clust_coeff(ConnMat);
disp(sprintf('Initial clustering coefficient (MIT toolbox) is: %f \n',C2));

% Degree distributions (Fig S3)
outgoing = sum(ConnMat,2)
incomming = sum(ConnMat,1)
HIST_O = histc(outgoing,1:1:ceil(max(outgoing)) ) ;
HIST_I = histc(incomming,1:1:ceil(max(incomming)) ) ;
bar(HIST_O);
figure;bar(HIST_I');


% Count neighbors: Probability distribution
for i=1:nPC
        kneighs(i) = length(kneighbors(ConnMat,i,1));
end
tempX = 1:1:ceil(max(kneighs));
HIST_KN = histc(kneighs,tempX ) ;
bar(HIST_KN);

% ConnMat = round(rand(nPC,nPC));





CC(1) = 0;
diffCC = 1;
t=2;
figure;
while 1 %diffCC > 0.001
% for t=1:1000
[ConnMat, probsMat] = reshapeNetwork(DistMat',ConnMat',NconnBins,NincomingProbs, NoutgoingProbs,recipBins,recipProbs); % problematic

% pause
% Histo = generateHisto(ConnMat', DistMat',allCombs', min(min(DistMat)), max(max(DistMat)));
% figure(2);bar(Histo);
% pause;
% close(2);
sortedCC = sort(clustCoeff(ConnMat),'descend');
cla;
plot(sortedCC);

CC(t) = nanmax(clustCoeff(ConnMat));
diffCC = abs(CC(t-1) - CC(t));
disp(sprintf('In time %d Average CC is %f, with diff %f\n',t,CC(t), CC(t-1) ) );
[C1,C2,C]=clust_coeff(ConnMat);
disp(sprintf('Clustering coefficient (MIT) is: %f \n',C2));
CC(t-1) = CC(t); 
t=t+1;

pause(0.1);
end

% Replicate figure 2 B, (Perin et al):
%Cluster of 3 pyramidals
for i=1:length(allCombs3)
    counter = 0;
    for j=1:length(idxCombs3)
        tempIdx = allCombs3(i,idxCombs3(j,:));
        counter = counter + ConnMat( tempIdx(1), tempIdx(2) );
    end
    cluster3(i) = counter;
end
clusterHist3 = histc(cluster3,0:length(idxCombs3)) ;
figure;plot(clusterHist3);figure(gcf);

%Cluster of 4 pyramidals
for i=1:length(allCombs4)
    counter = 0;
    for j=1:length(idxCombs4)
        tempIdx = allCombs4(i,idxCombs4(j,:));
        counter = counter + ConnMat( tempIdx(1), tempIdx(2) );
    end
    cluster4(i) = counter;
end
clusterHist4 = histc(cluster4,0:length(idxCombs4)) ;
figure;plot(clusterHist4);figure(gcf);

%Cluster of 5 pyramidals
for i=1:length(allCombs5)
    counter = 0;
    for j=1:length(idxCombs5)
        tempIdx = allCombs5(i,idxCombs5(j,:));
        counter = counter + ConnMat( tempIdx(1), tempIdx(2) );
    end
    cluster5(i) = counter;
end
clusterHist5 = histc(cluster5,0:length(idxCombs5)) ;
figure;plot(clusterHist5);figure(gcf);

%Cluster of 6 pyramidals
for i=1:length(allCombs6)
    counter = 0;
    for j=1:length(idxCombs6)
        tempIdx = allCombs6(i,idxCombs6(j,:));
        counter = counter + ConnMat( tempIdx(1), tempIdx(2) );
    end
    cluster6(i) = counter;
end
clusterHist6 = histc(cluster6,0:length(idxCombs6)) ;
figure;plot(clusterHist6);figure(gcf);

%%
% make diagonal ones (connected autapses for PCs):
ConnMat = ConnMat + eye(size(ConnMat)) ; 
ConnMat(find(ConnMat>1)) = 1 ;
% Pad connection matrix (for interneuron types):
ConnMat = padarray(ConnMat,[nAllCells-nPC, nAllCells-nPC],1,'post');
% Export parameter matrices in .hoc file:

% Write NMDA results to a .hoc file:
fid = fopen('importNetworkParameters.hoc','w');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');
fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref connMatrix, weightsMatrix\n');
fprintf(fid,'connMatrix = new Matrix(nAllCells, nAllCells)\n');
fprintf(fid,'weightsMatrix = new Matrix(nAllCells, nAllCells)\n');

fprintf(fid,'\n\n// Import parameters: (long-long text following!)\n\n');
% network connectivity:
for i=1:length(ConnMat)
    for j=1:length(ConnMat)
        fprintf(fid,'connMatrix.x[%d][%d] = %d\n',i-1,j-1, ConnMat(i,j));
    end
end
% Network synaptic weights
for i=1:length(ConnMat)
    for j=1:length(ConnMat)
        fprintf(fid,'weightsMatrix.x[%d][%d] = %f\n',i-1,j-1, ConnMat(i,j));
    end
end

fclose(fid);
%% OLD CODE
% figure;imagesc(ConnMat);

% Check single connectivity
TEMPMAT = (ConnMat .* DistMat);
TEMPMAT = TEMPMAT(TEMPMAT>0);
HIST = histc(TEMPMAT(:),connBins); % was connBins
SUMH = histc(DistMat(:),connBins);
HIST = HIST ./ SUMH;
% HIST = HIST / (length(TEMPMAT(:)));
bar(HIST);

clear TEMPMAT HIST SUMH
figure;
% Check unidirectional connectivity (reciprocal):
TEMPMAT = (ConnMat .* ConnMat');
TEMPMAT = (TEMPMAT .* DistMat);
TEMPMAT = TEMPMAT(TEMPMAT>0);
HIST = histc(TEMPMAT(:),recipBins); % was connBins
SUMH = histc(DistMat(:),recipBins);
HIST = HIST ./ SUMH;
bar(HIST);



% %Generate figure 3 C (Perin et al)
% for i=1:length(allCombs)
%     tempPerms = perms(allCombs(i,:));
%     tempPerms = unique(tempPerms(:,1:2),'rows');
%     for k=1:length(tempPerms)
%         tempDist(k) = DistMat( tempPerms(k,1), tempPerms(k,2) );
%     end
%     dist6(i) = mean(tempDist);
% end
% HIST6 = histc(dist6,[50,75,100,125,150,175,200,225]);
% HIST6 = HIST6 / (length(dist6));
