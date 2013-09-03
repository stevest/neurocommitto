close all;clear all;
Points = CreatePerinNetwork(28,15);
nPC = length(Points);
nPV = 0; % was 6
nCB = 0;
nCR = 0;
AllCells = [nPC , nPV , nCB , nCR];
nAllCells = sum(AllCells);

% DistMat = generateDistanceMat(Points(1:nPC,:)');
DistMat = generateDistanceMat(Points(1:nPC,:)',28*6,1);

% Use custom combinations function for effficiency:
% maxCombs = 10000;
allCombs3 = combntnsRND(nPC,3,4199);
allCombs4 = combntnsRND(nPC,4,8202);
allCombs5 = combntnsRND(nPC,5,11544);
allCombs6 = combntnsRND(nPC,6,12012);
idxCombs3 = [combntns(1:3,2); combntns(3:-1:1,2)];
idxCombs4 = [combntns(1:4,2); combntns(4:-1:1,2)];
idxCombs5 = [combntns(1:5,2); combntns(5:-1:1,2)];
idxCombs6 = [combntns(1:6,2); combntns(6:-1:1,2)];
%%
connBins = [20,50,90,125,160,190,225,260,295,340];
connProbs = [0.25, 0.21, 0.18, 0.15, 0.11, 0.09, 0.06, 0.04, 0.02, 0.01] ;
recipBins = [20,50,90,125,160,190,225,260,295,340];
recipProbs = [0.12, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.003]  ;

NconnBins = [0,1,2,3];
NincomingProbs = [0.12, 0.21, 0.25, 0.4] * 0.01;
NoutgoingProbs = [0.1, 0.25, 0.22, 0.18] * 0.01;

ConnMat = initializeNetwork(DistMat',connBins,connProbs,recipBins,recipProbs);
for i=1:1000
[ConnMat, probsMat] = reshapeNetwork(DistMat',ConnMat',NconnBins,NincomingProbs, NoutgoingProbs,recipBins,recipProbs); % problematic
end

% outgoing = sum(ConnMat,2)
% incomming = sum(ConnMat,1)
% HIST_O = histc(outgoing,0:1:ceil(max(outgoing)) ) ;
% HIST_I = histc(incomming,0:1:ceil(max(incomming)) ) ;
% bar(HIST_O);
% figure;bar(HIST_I');

% [ConnMat, probsMat] = reshapeNetwork(DistMat',ConnMat',NconnBins,NincomingProbs, NoutgoingProbs,recipBins,recipProbs); % problematic

% start=1;
% debug=[0,0,0];
% for i=1:2028
%     for j=start:2028
%         testCondition =  ConnMat(i,j) + ConnMat(j,i);
%         if(testCondition==0)
%             debug(1) = debug(1) + 1;
%         else if (testCondition==1)
%                 debug(2) = debug(2) + 1;
%             else
%                 debug(3) = debug(3) + 1;
%             end
%         end
%     end
%     start = start+1;
% end

start=1;
debug=[];
debugH = [];
ctr=1;
for i=1:length(ConnMat)
    for j=start:length(ConnMat)
        testCondition =  ConnMat(i,j) + ConnMat(j,i);
        
        if (testCondition==1)
                debug(ctr) = DistMat(i,j); 
                ctr=ctr+1;
        end
        
%         if (testCondition==2)
%                 debug(ctr) = DistMat(i,j); 
%                 ctr=ctr+1;
%         end
    end
    start = start+1;
end
debugH = histc(debug,[20,50,90,125,160,190,225,260,295,340] ) ;
figure;plot(debugH);figure(gcf);


% %Cluster of 6 pyramidals
% for i=1:length(allCombs6)
%     counter = 0;
%     for j=1:length(idxCombs6)
%         tempIdx = allCombs6(i,idxCombs6(j,:));
%         counter = counter + ConnMat( tempIdx(1), tempIdx(2) );
%     end
%     cluster6(i) = counter;
% end
% clusterHist6 = histc(cluster6,0:length(idxCombs6)) ;
% figure;plot(clusterHist6);figure(gcf);
