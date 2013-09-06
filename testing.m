close all;clear all;

mex generateDistanceMat.c
mex initializeNetwork.c
mex reshapeNetwork.c

% Initialize somata positions as in Perin et al.
Points = CreatePerinNetwork(28,15);

% Generate distance matrix:
DistMat = generateDistanceMat(Points(1:length(Points),:)', 1, 28*12);

%% Main

% Import initial connectivity probabilities from Perin et al. (figure 1)
connBins = [20,50,90,125,160,190,225,260,295,340];
connProbs = [0.25, 0.21, 0.18, 0.15, 0.11, 0.09, 0.06, 0.04, 0.02, 0.01];
recipBins = [20,50,90,125,160,190,225,260,295,340];
recipProbs = [0.12, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.003];

% Update probabilities from Perin et al. (Figure S4)
NconnBins = [0,1,2,3];
NincomingProbs = [0, 0.33, 0.66, 0.99] ;
NoutgoingProbs = [0, 0.33, 0.66, 0.99] ;

% Initialize network: Create connectivity based on intersomatic distance: 
ConnMat = initializeNetwork(DistMat',connBins,connProbs,recipBins,recipProbs);

% Iterative step
for t=1:100
[ConnMat, probsMat] = reshapeNetwork(DistMat',ConnMat',NconnBins,NincomingProbs, NoutgoingProbs,recipBins,recipProbs); % Probabilities map(probsMat) is optional
[~,~,CC(t)] = clust_coeff(ConnMat) ;
end


%% Following code to recreate functions:

% Use custom combinations function for effficiency:
allCombs3 = combntnsRND(length(Points),3,4199);
allCombs4 = combntnsRND(length(Points),4,8202);
allCombs5 = combntnsRND(length(Points),5,11544);
allCombs6 = combntnsRND(length(Points),6,12012);
idxCombs3 = [combntns(1:3,2); combntns(3:-1:1,2)];
idxCombs4 = [combntns(1:4,2); combntns(4:-1:1,2)];
idxCombs5 = [combntns(1:5,2); combntns(5:-1:1,2)];
idxCombs6 = [combntns(1:6,2); combntns(6:-1:1,2)];

% Code to recreate Fig S3: (fail)
outgoing = sum(ConnMat,2);
incomming = sum(ConnMat,1);
HIST_O = histc(outgoing,0:1:ceil(max(outgoing)) ) ;
HIST_I = histc(incomming,0:1:ceil(max(incomming)) ) ;
bar(HIST_O);
figure;bar(HIST_I');


% Code to recreate Fig 5E: (fail)
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
clusterHist3 = clusterHist3/norm(clusterHist3);
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
