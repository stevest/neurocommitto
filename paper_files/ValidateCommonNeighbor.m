close all;clear all; clc;
aa=1;
cb = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.7773    0.9102    1.0000;
    1.0000    0.8555    0.7969;
    1.0000    0.8867    0.6289];
rng('default')
rng(aa)
total_prob = @(x) 0.22.*exp(-0.0052*x);
recipProbsPC2PC = @(x) 0.12 .* exp(-0.0085* x);
connProbsPC2PC =@(x) ((0.22.*exp(-0.0052*x)) - (0.12 .* exp(-0.0085* x)))*2;

% Number of neurons
nPC = 75%300;
cubeRange = 200;%100:100:1000;
many_nets=10;
PC2PC_str = zeros(nPC,nPC,many_nets);
PC2PC_d   = zeros(nPC,nPC,many_nets);
PC2PC_rnd = zeros(nPC,nPC,many_nets);
distPC2PC = ones(nPC,nPC,many_nets);

for cdim = 1:length(cubeRange) % find(cubeRange==200) %
    cube_dimensions = cubeRange(cdim);
    for k=1:many_nets
        k
        % Change network in space and check what is happening:
        [PCsomata, distPC2PC(:,:,k,cdim)]= CreateRandomNetwork(nPC, cube_dimensions);
        [PC2PC_d(:,:,k,cdim), PC2PC_str(:,:,k,cdim), ~, ~, ~, ~]= create_graph_CN(distPC2PC(:,:,k,cdim),connProbsPC2PC,recipProbsPC2PC, total_prob);
        prob_conn_rnd_ind =0.13;
        [PC2PC_rnd(:,:,k,cdim),~, ~, ~] = create_graph_WS(nPC,prob_conn_rnd_ind,'ind'); % 1.0=Random graph, 0.0=Watts-Strogatz graph
    end
    % Check common neighbour if sufficient assumption:
end

%%

averageConns_str = cell(6,length(cubeRange));
averageConns_rnd = cell(6,length(cubeRange));
averageConns_d = cell(6,length(cubeRange));
clusterSize = [3,4,5,6,7,8];
for cdim = 1:length(cubeRange)
    for k=1:length(clusterSize)
        [cdim k]
        for randomClusters=1:50000
            PCpermut = randperm(nPC,clusterSize(k));
            PCbatch_str = PC2PC_str(PCpermut,PCpermut,ceil(rand(1)*many_nets),cdim);
            PCbatch_rnd = PC2PC_rnd(PCpermut,PCpermut,ceil(rand(1)*many_nets),cdim);
            PCbatch_d = PC2PC_d(PCpermut,PCpermut,ceil(rand(1)*many_nets),cdim);
            averageConns_str{k,cdim} = [averageConns_str{k,cdim}, sum(PCbatch_str(:))];
            averageConns_rnd{k,cdim} = [averageConns_rnd{k,cdim}, sum(PCbatch_rnd(:))];
            averageConns_d{k,cdim} = [averageConns_d{k,cdim}, sum(PCbatch_d(:))];
        end
    end
end

%%

connsRange = [0:1:64];
z_alpha = 1.96;
h = figure;hold on;
k=6; % Cluster size
for q = 1
    connsProb_str = histc(averageConns_str{k,q},connsRange);
    connsProb_str = connsProb_str/randomClusters;
    plot(connsProb_str,'Color',[ 0.8500    0.3250    0.0980 ],'Linewidth',2);
    connsProb_rnd = histc(averageConns_rnd{k,q},connsRange);
    connsProb_rnd = connsProb_rnd/randomClusters;
    plot(connsProb_rnd,'Color',cb(1,:),'Linewidth',2);
    connsProb_d = histc(averageConns_d{k,q},connsRange);
    connsProb_d = connsProb_d/randomClusters;
    plot(connsProb_d,'Color',[0,0,0],'Linewidth',2);
    n = randomClusters;
    significance = zeros(1,length(connsProb_str));
    % Kolmogorov-Smirnov does not applay to discrete distributions:
% http://stats.stackexchange.com/questions/1047/is-kolmogorov-smirnov-test-valid-with-discrete-distributions
z = zeros(1,length(connsProb_str));
%     for k = 1:length(connsProb_str)
%         p1 = connsProb_str(k) ;
%         p2 = connsProb_d(k) ;
%         p = (n*p1 + n*p2) / (n+n) ;
%         z(k) = (p1-p2) / sqrt(p*(1-p)*(2/n)) ; 
%         x1 = averageConns_str{q} == k;
%         x2 = averageConns_d{q} == k;
%         if abs(z(k)) > z_alpha
%             significance(k) = 1;
%             plot(k,connsProb_str(k)+connsProb_str(k),'k*');
%         end
%     end
    
    set(gca,'YScale','log');
    xlabel(sprintf('Connections in clusters of %d cells.',clusterSize(k)));
    ylabel('Frequency');
    legend({'Structured','Random','Expected (distance)'});
%     pause();
%     cla;
end
% [h,p] = kstest2(averageConns_str{q},averageConns_d{q})



