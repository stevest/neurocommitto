close all;clear all;
cd('C:\Users\steve\Documents\GitHub\neurocommitto');

% mex generateDistanceMat.c
% mex initializeNetwork.c
% mex  reshapeNetwork.c
% mex generateHisto.c
% mex -g convergeConnectivity.c

% Initialize somata positions as in Perin et al.
Points = CreatePerinNetwork(28,15, 2);
scatter3(Points(1:144,1),Points(1:144,2), Points(1:144,3));
axis equal;
% Generate distance matrix:
DistMat = generateDistanceMat(Points', 0, 28*12);

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
% Histo of distances:
[a,b] = hist(DistMat(:), connBins);
a=a/max(a);
ConnMat = initializeNetwork(DistMat',a,connBins,connProbs,recipBins,recipProbs);
totalInput = 10;

tic;
[ConnMat, probsMat] = reshapeNetwork(DistMat',ConnMat',a,NconnBins,NincomingProbs, NoutgoingProbs,recipBins,recipProbs, totalInput); % Probabilities map(probsMat) is optional
toc
imagesc(ConnMat)

Layers = 1:144:length(Points);
layerLen = 143;

% for only 12 cells of given layer:
dozen = floor((rand(144,1) * 143)+Layers(ceil(rand(1)*length(Layers)))) ;

% Use custom combinations function for effficiency:
allCombs3 = combntns(1:12,3);
allCombs4 = combntns(1:12,4);
allCombs5 = combntns(1:12,5);
allCombs6 = combntns(1:12,6);
allCombs7 = combntns(1:12,7);
allCombs8 = combntns(1:12,8);
idxCombs3 = [combntns(1:3,2); combntns(3:-1:1,2)];
idxCombs4 = [combntns(1:4,2); combntns(4:-1:1,2)];
idxCombs5 = [combntns(1:5,2); combntns(5:-1:1,2)];
idxCombs6 = [combntns(1:6,2); combntns(6:-1:1,2)];
idxCombs7 = [combntns(1:7,2); combntns(7:-1:1,2)];
idxCombs8 = [combntns(1:8,2); combntns(8:-1:1,2)];

% Code to recreate Fig 5E: (fail)
%Cluster of 3 pyramidals
for t = 1:1000
    dozen = floor((rand(144,1) * 143)+Layers(ceil(rand(1)*length(Layers)))) ;
    for i=1:length(allCombs3)
        counter = 0;
        for j=1:length(idxCombs3)
            tempIdx = allCombs3(i,idxCombs3(j,:));
            counter = counter + ConnMat( dozen(tempIdx(1)), dozen(tempIdx(2)) );
        end
        cluster3(t,i) = counter;
    end
    clusterHist3(t,:) = histc(cluster3(t,:),0:length(idxCombs3)) ;
    clusterHist3(t,:) = clusterHist3(t,:)/norm(clusterHist3(t,:));
end
[aa,b] = hist(DistMat(:), length(idxCombs3)+1);
aa=aa/max(aa);
figure;plot(mean(clusterHist3)./aa);figure(gcf);


%Cluster of 4 pyramidals
for t = 1:1000
    dozen = floor((rand(144,1) * 143)+Layers(ceil(rand(1)*length(Layers)))) ;
    for i=1:length(allCombs4)
        counter = 0;
        for j=1:length(idxCombs4)
            tempIdx = allCombs4(i,idxCombs4(j,:));
            counter = counter + ConnMat( dozen(tempIdx(1)), dozen(tempIdx(2)) );
        end
        cluster4(t,i) = counter;
    end
    clusterHist4(t,:) = histc(cluster4(t,:),0:length(idxCombs4)) ;
    clusterHist4(t,:) = clusterHist4(t,:)/norm(clusterHist4(t,:));
end
[aa,b] = hist(DistMat(:), length(idxCombs4)+1);
aa=aa/max(aa);
figure;plot(mean(clusterHist4));figure(gcf);



%Cluster of 6 pyramidals
for t = 1:1000
    dozen = floor((rand(144,1) * 143)+Layers(ceil(rand(1)*length(Layers)))) ;
    for i=1:length(allCombs6)
        counter = 0;
        for j=1:length(idxCombs6)
            tempIdx = allCombs6(i,idxCombs6(j,:));
            counter = counter + ConnMat( dozen(tempIdx(1)), dozen(tempIdx(2)) );
        end
        cluster6(t,i) = counter;
    end
    clusterHist6(t,:) = histc(cluster6(t,:),0:length(idxCombs6)) ;
    clusterHist6(t,:) = clusterHist6(t,:)/norm(clusterHist6(t,:));
end
[aa,b] = hist(DistMat(:), length(idxCombs6)+1);
aa=aa/max(aa);
figure;plot(mean(clusterHist6)./aa);figure(gcf);


%Motifs of 6 pyramidals
mot6=[];
motBins = [50, 75, 100, 125, 150, 175, 200, 225];
for m = 1:length(motBins)
    motif6 = zeros(1000,length(allCombs6) );
    for t = 1:100
        dozen = floor((rand(144,1) * 143)+Layers(ceil(rand(1)*length(Layers)))) ;
        for i=1:length(allCombs6)
            avgDist = 0;
            for j=1:length(idxCombs6)
                tempIdx = allCombs6(i,idxCombs6(j,:));
                avgDist = avgDist + DistMat( dozen(tempIdx(1)), dozen(tempIdx(2)) );
            end
            avgDist = avgDist / length(idxCombs6);
            if((avgDist>(motBins(m)-25)) && (avgDist<motBins(m)))
                counter=0;
                for j=1:length(idxCombs6)
                    tempIdx = allCombs6(i,idxCombs6(j,:));
                    counter = counter + ConnMat( dozen(tempIdx(1)), dozen(tempIdx(2)) );
                end
                motif6(t,i) = counter;
            end
        end
    end
    
    mot6(m) = mean(motif6(find(motif6>0)));
end

%%


dta = []
bdta = []
for i=1:length(DistMat)
    for j=1:i
        if (ConnMat(i,j)>0)
            dta(end+1) = DistMat(i,j);
            if (ConnMat(j,i)>0)
                bdta(end+1) = DistMat(j,i);
            end
        end
    end
end
% figure()
hist( dta(:) , 10 )

% [CS,CR] = generateHisto(ConnMat.*DistMat,connBins)
%
% % CSp = CS/norm(CS);
% CSp = CS*(1/sum(CS));
% % CRp = CR/norm(CR);
% CRp = CRp*(1/sum(CRp));
%
% CS = CSp./CS;
% CR = CRp./CR;
%
%
% bar(CR);
% bar(CS);

%% Initialize network: Create connectivity based on intersomatic distance:
% ConnMat = initializeNetwork(DistMat',connBins,connProbs,recipBins,recipProbs);
ConnMat = ones(2028,2028);
pntsPerTime = 20;

[CS,CR] = generateHisto(ConnMat.*DistMat,connBins);
CS = CS/norm(CS);
CR = CR/norm(CR);
cvgS = abs(sum(connProbs - CS));
cvgR = abs(sum(recipProbs - CR));

srand = clock;
rng(srand(end));
ctr = 1;
for t = 1:200000
    ri = ceil(rand(1,pntsPerTime)*2028);
    rj = ceil(rand(1,pntsPerTime)*2028);
    
    % flipedEls = xor(diag(ConnMat(ri,rj)),ones(pntsPerTime,1));
    for i=1:pntsPerTime
        ConnMat(ri(i),rj(i)) = ~ConnMat(ri(i),rj(i));
    end
    
    [CS,CR] = generateHisto(ConnMat.*DistMat,connBins);
    CS = CS/norm(CS);
    CR = CR/norm(CR);
    
    if( (abs(sum(connProbs - CS)) > cvgS) || (abs(sum(recipProbs - CR)) > cvgR) )
        for i=1:pntsPerTime
            ConnMat(ri(i),rj(i)) = ~ConnMat(ri(i),rj(i));
        end
    else
        improvementS(ctr) = cvgS - abs(sum(connProbs - CS));
        improvementR(ctr) = cvgR - abs(sum(recipProbs - CR));
        disp(sprintf('@%d Improvement=%f,ctr=%d error=%f\n',t,mean([improvementS(ctr),improvementR(ctr)]),ctr,mean([cvgS,cvgR])) );
        ctr = ctr +1;
        cvgS = abs(sum(connProbs - CS));
        cvgR = abs(sum(recipProbs - CR));
    end
    
    
end

%% Iterative step
for t=1:100
    [ConnMat, probsMat] = reshapeNetwork(DistMat',ConnMat',a,NconnBins,NincomingProbs, NoutgoingProbs,recipBins,recipProbs); % Probabilities map(probsMat) is optional
    [~,~,CC(t)] = clust_coeff(ConnMat) ;
end


%% Following code to recreate functions:
% % Stirling's Approximation
% n = 20;
% k = 8;
% d = n-k;
%
% Stirling = ceil(((n^n)*exp(-n)*sqrt(2*pi*n)) / (((k^k)*exp(-k)*sqrt(2*pi*k)) * ((d^d)*exp(-d)*sqrt(2*pi*d))));
% CBN = cc_combntns(Stirling+100);

% Use custom combinations function for effficiency:
allCombs3 = combntns(Layers(1):Layers(1)+layerLen,3);
allCombs4 = combntns(Layers(1):Layers(1)+layerLen,4);
allCombs5 = combntns(Layers(1):Layers(1)+layerLen,5);
allCombs6 = combntns(Layers(1):Layers(1)+layerLen,6);
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
clusterHist5 = clusterHist5/norm(clusterHist5);
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
