close all;
snormal = @(x,m,a,sigma) ((1+erf((a.*x)/sqrt(2)))/2) .* normpdf(x,m,sigma);
rang = 0:0.01:8;
maxNeighbors = 25;
alpha = 3*ones(1,maxNeighbors);%linspace(3,0.1,maxNeighbors)%.*exp(linspace(0,0.5,maxNeighbors)).^-2;%[20,10,2,1];
maxSigma = 5 ;
minSigma = 0.1 ;
expClimb = (exp(linspace(0,0.5,maxNeighbors)).^5)-1;
% figure;plot(expClimb);
expClimb = (expClimb/max(expClimb)*maxSigma)+minSigma;
m = ones(1,maxNeighbors).*expClimb.*2;%exp(linspace(0,0.6,maxNeighbors))-1;
m = m-m(1);
sigma = 0.5*ones(1,maxNeighbors).*expClimb%.*(1./(1+exp(-linspace(-2,5,maxNeighbors))))%.*exp(linspace(0,0.5,maxNeighbors)).^-10;%[0.15, 0.2, 4, 6];
normFactors = ones(1,maxNeighbors);%[0.0045, 0.006, 0.02, 0.02];
% LogNormal PDF parameters:
parmhat(1,1:maxNeighbors) = linspace(0.1,10,maxNeighbors);
parmhat(2,1:maxNeighbors) = 2.6;
CNsnormal = [];
% Probabilities changed to replicate mean EPSP amplitude per connections in
% cluster (Perin et al, Fig.6A):
maxCDF = zeros(maxNeighbors,length(rang));
for k=1:maxNeighbors
    CNsnormal(k,:) = snormal(rang, m(k), alpha(k), sigma(k))*normFactors(k);
    maxCDF(k,:) = cumsum(CNsnormal(k,:))/max(cumsum(CNsnormal(k,:)));
%     CNsnormal(k,:) = lognpdf(rang,parmhat(1,k),parmhat(2,k));
end
figure;
plot(CNsnormal(:,:)');hold on;
set(gca,'xtick',1:60:length(rang));
set(gca,'xticklabel',rang(1:60:end));
%%
% many_nets = 100;
% maxCDF = [max(cumsum(CNsnormal(1,:))), max(cumsum(CNsnormal(2,:))), max(cumsum(CNsnormal(3,:))),  max(cumsum(CNsnormal(4,:)))];
% pre compute the weights per cn distribution:
maxPrecomps = 75*75*4;
sterngthArr = zeros(maxNeighbors,maxPrecomps);
for n = 1:maxNeighbors
    n
    ctr = 1;
    while ctr <= maxPrecomps%length(find(CN_str==k-1))
        randomSampling = rand(1);
%         strength = rang(find(histc(randomSampling * max(cumsum(CNsnormal(k,:))),cumsum(CNsnormal(k,:)))));
        strength = rang(find(histc(randomSampling,maxCDF(n,:))));
        if ~isempty(strength)
            sterngthArr(n,ctr) = strength(1);
            ctr = ctr+1;
        end
    end
end
%%
figure;hold on;
for n = 1:maxNeighbors
    n
    plot(histc(sterngthArr(n,:),rang))
end
set(gca,'xtick',1:60:length(rang));
set(gca,'xticklabel',rang(1:60:end));    

weights_str = ones(75,75,many_nets);
for stateID = 1:many_nets
    stateID
    tmp_weights_str = ones(75,75);
    [CN_str] = m_commonNeighbors(PC2PC_str(1:75,1:75,stateID));
    CN_str(CN_str>maxNeighbors) = maxNeighbors;
    for k=1:maxNeighbors
%         sterngthArr = [];
%         ctr = 0;
%         while ctr < numel(CN_str)*2%length(find(CN_str==k-1))
%             randomSampling = rand(1);
% %             strength = rang(find(histc(rand(1)*max(cumsum(CNsnormal(k,:))),cumsum(CNsnormal(k,:)))));
%             strength = rang(find(histc(randomSampling * max(cumsum(CNsnormal(k,:))),cumsum(CNsnormal(k,:)))));
%             if ~isempty(strength)
%                 sterngthArr = [sterngthArr, strength(1)];
%                 ctr = ctr+1;
%             end
%         end
        tmp_weights_str(CN_str==k-1) = sterngthArr(k,randperm(maxPrecomps,sum(sum(CN_str==k-1))));
    end
    weights_str(:,:,stateID) = tmp_weights_str .* PC2PC_str(1:75,1:75,stateID);
end


averageEPSP = cell(many_nets,30);
for kk = 1:many_nets
    kk
    for k=1:1000
        PCpermut = randperm(75,6);
        tmpNetId = 1;%ceil(rand(1)*many_nets);
        PCbatch = PC2PC_str(PCpermut,PCpermut,tmpNetId);
        EPSPs = weights_str(PCpermut,PCpermut,tmpNetId);
        EPSPs(EPSPs==0) = [];
        if ~isempty(EPSPs)
            averageEPSP{kk,sum(PCbatch(:))+1} = [averageEPSP{kk,sum(PCbatch(:))+1}, EPSPs];
        end
    end
end
AvrEPSPamplitude = cell(many_nets,1);
for kk = 1:many_nets
    kk
    AvrEPSPamplitude{kk,1} = cell2mat(cellfun(@(x) mean(x),averageEPSP(kk,:),'Uniformoutput',false));
end
EPSPdata = cell2mat(AvrEPSPamplitude);
EPSPstd = nanstd(EPSPdata,1);
EPSPmean = nanmean(EPSPdata,1);
EPSPsem = EPSPstd./sqrt( sum(~isnan(EPSPdata),1) );
figure;errorbar(EPSPmean,EPSPsem);


averageEPSPhist = cellfun(@(x) histc(x,[0:0.1:6]),averageEPSP(1,:),'UniformOutput', false);
figure;
for k=1:length(averageEPSPhist)
    subplot(1,length(averageEPSPhist),k);box off;
    barh(averageEPSPhist{k});set(gca,'xticklabel','','yticklabel','');
end

%% validate cluster weights:
tmpWeightBins = cell(1,maxNeighbors);
for k=1:maxNeighbors
        tmpWeightBins{k} = zeros(1,length(rang));
end
for stateID = 1:many_nets
    stateID
    [CN_str] = m_commonNeighbors(PC2PC_str(1:75,1:75,stateID));
    tmpWeights = weights_str(:,:,stateID);
    for k=1:maxNeighbors
        tmpWeightBins{k} = tmpWeightBins{k} + histc(tmpWeights(CN_str==k),rang)';
    end
end

plot(cell2mat(tmpWeightBins')')


%% Find clusters in networks
WB_str = {};
cellsPerCluster_str = cell(1,many_nets);
cellsPerCluster_rnd = cell(1,many_nets);
% WB_str_BC = WB_str;
WB_rnd = {};
% WB_rnd_BC = WB_rnd;
for stateID = 1:many_nets
    % Find nearest neighbors.
    [CN_str] = m_commonNeighbors(PC2PC_str(:,:,stateID));
%     [ CN_d ] = m_commonNeighbors(PC2PC_d(:,:,stateID));
    [CN_rnd] = m_commonNeighbors(PC2PC_rnd(:,:,stateID));
    
    % keep only the pairs that are connected.
    mergedCN_str = CN_str .* PC2PC_str(:,:,stateID);
%     mergedCN_d   =  CN_d  .* PC2PC_d(:,:,stateID);
    mergedCN_rnd = CN_rnd .* PC2PC_rnd(:,:,stateID);
    
    % Affinity propagation algorithm (Frey, Dueck, 2007):
    % Force different # of cluster to get as many clusters as you
    % can from the populations with high min cells per cluster:
    ClNo = 5;
    minClustSize = 9;
    Flg = 1;
    while Flg
%         for i=1:5
            labels_str=[];
            NC_str=[];
            % Nassi, Changed in apclusterK total # of bisections from 20 to 50 (line 75)
            [idx_str,~,~,~,pref_str]=apclusterK(mergedCN_str, ClNo,0)
            [targ_str,~,labels_str] = unique(idx_str);
            if ( min(histc(labels_str,1:ClNo)') >= minClustSize)
                Flg = 0;
                ClNo
                break;
            end
%         end
%         ClNo = ClNo - 1;
    end
    cl_labels_str(:,stateID)=labels_str;
    [within_str(stateID),between_str(stateID),WB_str(stateID)] = calculateWithinBetween(PC2PC_str(:,:,stateID),labels_str);
%     [within_w_str(stateID),between_w_str(stateID),WB_w_str(stateID)] = calculateWithinBetween(weights_str(:,:,stateID),labels_str);
    cellsPerCluster_str(stateID) = {histc(labels_str,1:size(targ_str,1))'};

    Flg = 1;
    while Flg
%         for i=1:5
            labels_rnd=[];
            NC_rnd=[];
            % Nassi, Changed in apclusterK total # of bisections from 20 to 50 (line 75)
            [idx_rnd,~,~,~,pref_rnd]=apclusterK(mergedCN_rnd, ClNo,0)
            [targ_rnd,~,labels_rnd] = unique(idx_rnd);
            if ( min(histc(labels_rnd,1:ClNo)') >= minClustSize)
                Flg = 0;
                ClNo
                break;
            end
%         end
%         ClNo = ClNo - 1;
    end
    cl_labels_rnd(:,stateID)=labels_rnd;
    [within_rnd(stateID),between_rnd(stateID),WB_rnd(stateID)] = calculateWithinBetween(PC2PC_rnd(:,:,stateID),labels_rnd);
%     [within_w_str(stateID),between_w_str(stateID),WB_w_str(stateID)] = calculateWithinBetween(weights_str(:,:,stateID),labels_str);
    cellsPerCluster_rnd(stateID) = {histc(labels_rnd,1:size(targ_rnd,1))'};
    
end
% Check weights distribution per cluster:
averageClusterEPSP = cell(many_nets,10);
for stateID=1:many_nets
    for cl = 1:length(cellsPerCluster_str{stateID})
        idx = find(cl_labels_str(:,1)==cl);
        PCbatch = PC2PC_str(idx,idx,stateID);
        EPSPs = weights_str(idx,idx,stateID);
        EPSPs(EPSPs==0) = [];
        if ~isempty(EPSPs)
            averageClusterEPSP{stateID,cl} = [averageClusterEPSP{stateID,cl}, EPSPs];
        end
    end
end
rangmin = rang(1:10:end)
histoClusterEPSP = cellfun(@(x) histc(x,rangmin),averageClusterEPSP,'uniformoutput',false);
histoClusterEPSP = reshape(histoClusterEPSP,1,[]);
histoClusterEPSP(cellfun(@isempty,histoClusterEPSP)) = [];
EPSPdata = cell2mat(histoClusterEPSP');
EPSPstd = nanstd(EPSPdata,1);
EPSPmean = nanmean(EPSPdata,1);
EPSPsem = EPSPstd./sqrt( sum(~isnan(EPSPdata),1) );
figure;errorbar(EPSPmean,EPSPsem);
set(gca,'xtick',1:60:length(rangmin));
set(gca,'xticklabel',rangmin(1:60:end));  
figure;
h=plot(EPSPdata');
set(h,'alpha',0.2)

set(gca,'xtick',1:60:length(rangmin));
set(gca,'xticklabel',rangmin(1:60:end));  

%%
snormal = @(x,m,a,sigma) ((1+erf((a.*x)/sqrt(2)))/2) .* normpdf(x,m,sigma);
rang = 0:0.01:8;
m = [0,0,0,0];
alpha = [20,10,5,5];
sigma = [0.45, 0.55, 1.3, 2];
CNsnormal = [];
CNsnormal(1,:) = snormal(rang, m(1), alpha(1), sigma(1))*0.012;
CNsnormal(2,:) = snormal(rang, m(2), alpha(2), sigma(2))*0.014;
CNsnormal(3,:) = snormal(rang, m(3), alpha(3), sigma(3))*0.02;
CNsnormal(4,:) = snormal(rang, m(4), alpha(4), sigma(4))*0.02;
% data = snormal(rang, m(1), alpha(1), sigma(1))*0.012;
% [parmhat,parmci] = lognfit(data);
% plot(data);hold on;
% plot(lognpdf(rang,parmhat(1),parmhat(2)));
% CNsnormal(1,:) = lognpdf(rang,parmhat(1),parmhat(2));
% [parmhat,parmci] = lognfit(snormal(rang, m(2), alpha(2), sigma(2))*0.014);
% CNsnormal(2,:) = lognpdf(rang,parmhat(1),parmhat(2));
% [parmhat,parmci] = lognfit(snormal(rang, m(3), alpha(3), sigma(3))*0.02);
% CNsnormal(3,:) = lognpdf(rang,parmhat(1),parmhat(2));
% [parmhat,parmci] = lognfit(snormal(rang, m(4), alpha(4), sigma(4))*0.02);
% CNsnormal(4,:) = lognpdf(rang,parmhat(1),parmhat(2));
figure;
plot(CNsnormal(1,:),'b');hold on;
plot(CNsnormal(2,:),'r');hold on;
plot(CNsnormal(3,:),'g');hold on;
plot(CNsnormal(4,:),'k');hold on;
set(gca,'xtick',1:60:length(rang));
set(gca,'xticklabel',rang(1:60:end));

%%
sterngthArr = cell(1,4);
maxCDF = [0.5815, 0.6624, 0.9529,  0.9692];
% maxCDF = [3, 3, 1, 1];
CDF = [cumsum(CNsnormal(1,:));cumsum(CNsnormal(2,:));cumsum(CNsnormal(3,:));cumsum(CNsnormal(4,:))];
% for k=1:4
%     CDF(k,:) = CDF(k,:)/max(CDF(k,:));
% end
% plot(CDF')
for k=1:4
    ctr = 1;
    while ctr < 10000%length(find(CN_str==k-1))
        randomSampling = rand(1);
        if randomSampling > maxCDF(k)
            ctr = ctr+1
            continue;
        end
        strength = rang(find(histc(randomSampling,CDF(k,:))));
        if ~isempty(strength)
            sterngthArr{k} = [sterngthArr{k}, strength(1)];
            ctr = ctr+1
        end
    end
end

%%
factor = 0.0001
figure;hold on;
plot(CNsnormal(1,:),'b');hold on;
tmp = histc(sterngthArr{1},rang);
plot(tmp*factor,'b');
plot(CNsnormal(2,:),'r');hold on;
tmp = histc(sterngthArr{2},rang);
plot(tmp*factor,'r');
plot(CNsnormal(3,:),'g');hold on;
tmp = histc(sterngthArr{3},rang);
plot(tmp*factor,'g');
plot(CNsnormal(4,:),'k');hold on;
tmp = histc(sterngthArr{4},rang);
plot(tmp*factor,'k');
set(gca,'xtick',1:60:length(rang));
set(gca,'xticklabel',rang(1:60:end));

%%
factor = 0.0001
figure;hold on;
plot(cumsum(CNsnormal(1,:)),'b');hold on;
tmp = histc(sterngthArr{1},rang);
plot(cumsum(tmp*factor),'b');
plot(cumsum(CNsnormal(2,:)),'r');hold on;
tmp = histc(sterngthArr{2},rang);
plot(cumsum(tmp*factor),'r');
plot(cumsum(CNsnormal(3,:)),'g');hold on;
tmp = histc(sterngthArr{3},rang);
plot(cumsum(tmp*factor),'g');
plot(cumsum(CNsnormal(4,:)),'k');hold on;
tmp = histc(sterngthArr{4},rang);
plot(cumsum(tmp*factor),'k');
set(gca,'xtick',1:60:length(rang));
set(gca,'xticklabel',rang(1:60:end));


