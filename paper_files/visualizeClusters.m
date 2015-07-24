%% Distance dependent Network of different sizes:
total_prob = @(x) 0.22.*exp(-0.0052*x);
recipProbsPC2PC = @(x) 0.12 .* exp(-0.0085* x);
connProbsPC2PC =@(x) ((0.22.*exp(-0.0052*x)) - (0.12 .* exp(-0.0085* x)))*2;

nCells = 100;
cube_dimensions=200;
density = nCells/(cube_dimensions^3);

brainbox_size = 200:100:1000;
cellsPerCluster_d = {};
cm = lines(100);
for k=1:length(brainbox_size)
    nCells = density * (brainbox_size(k)^3);
    % 3-d points in space. Arguments: #of cells, max seperation distance.
    [PCsomata, distPC2PC]= CreateRandomNetwork(nCells, brainbox_size(k));
    [PC2PC_d, ~, ~, ~, ~, pp(2:3)]= create_graph_CN(distPC2PC,connProbsPC2PC,recipProbsPC2PC, total_prob);
    
    % Find nearest neighbors.
    [CN_d] = m_commonNeighbors(PC2PC_d);
    
    % keep only the pairs that are connected.
    mergedCN_d = CN_d .* PC2PC_d;

    ClNo = 5;
    while 1
        labels_d=[];
        [idx_d,~,~,~,~]=apclusterK(mergedCN_d, ClNo,0)
% [idx_str,~,~,~,pref_str]=apcluster(mergedCN_str, ClNo,'dampfact',0.9, ...
%             'convits',200,'maxits',2000,'nonoise');
        [targ_d,~,labels_d] = unique(idx_d);
        if ( min(histc(labels_d,1:ClNo)') >= minClustSize)
            break;
        end
    end
    cl_labels_d(:,stateID)=labels_d;
%     [within_str(stateID),between_str(stateID),WB_str(stateID)] = calculateWithinBetween(PC2PC_str(:,:,stateID),labels_str);
    NC_d = length(unique(labels_d));
    cellsPerCluster_d(k) = {histc(labels_d,1:size(targ_d,1))'};
    % Plot the result:
    figure;hold on;
    for c = 1:NC_d
        cells = find(labels_d == c);
        scatter3(PCsomata(cells,1),PCsomata(cells,2),PCsomata(cells,3),60,cm(c,:),'filled');
    end

end



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

many_nets=100;

PC2PC_str = zeros(nPC,nPC,many_nets);
PC2PC_d   = zeros(nPC,nPC,many_nets);
PC2PC_rnd = zeros(nPC,nPC,many_nets);

distPC2PC = ones(75,75,many_nets);
for k=1:many_nets
    k
    % Change network in space and check what is happening:
    [PCsomata, distPC2PC(:,:,k)]= CreateRandomNetwork(nPC, cube_dimensions);
    [PC2PC_d(:,:,k), PC2PC_str(:,:,k), coeffs_global(k, 2:3), coeffs_local(k,2:3), ~, pp(k,2:3)]= create_graph_CN(distPC2PC(:,:,k),connProbsPC2PC,recipProbsPC2PC, total_prob);
    prob_conn_rnd_ind =pp(1,3);%0.13; 
    % na to tre3w polles fores k na kratisw afto pou einai pio konta sto
    % input independent probability...
    [PC2PC_rnd(:,:,k), coeffs_global(k, 1), coeffs_local(k,1), pp(k, 1)] = create_graph_WS(nPC,prob_conn_rnd_ind,'ind'); % 1.0=Random graph, 0.0=Watts-Strogatz graph
end


%% Visualise big network clusters:
cm = lines(run.NC_str);
clusters = randperm(run.NC_str,10);
figure;hold on;
for c = 1:run.NC_str
    cells = find(run.labels_str == c);
    scatter3(PCsomata(cells,1),PCsomata(cells,2),PCsomata(cells,3),60,[0.7,0.7,0.7],'filled');
end
for c = clusters;%1:run.NC_rnd
    cells = find(run.labels_str == c);
    scatter3(PCsomata(cells,1),PCsomata(cells,2),PCsomata(cells,3),60,cm(c,:),'filled');
    tmp = run.state_str(cells,cells);
    for k = 1:length(cells)
        for j = 1:length(cells)
            if k ~= j
               plot3([PCsomata(cells(k),1),PCsomata(cells(j),1)],[PCsomata(cells(k),2),PCsomata(cells(j),2)],[PCsomata(cells(k),3),PCsomata(cells(j),3)],'Color',cm(c,:)); 
            end
        end
    end
end

axis equal;


%% PCA plot
X = CN_d;
X(X>3) = 3;
X = X - max(X(:));
D = pdist(X,'euclidean');
[Y,e] = cmdscale(D);

figure;hold on;plot(e);plot(e,'*');
xlabel('MDS Eigenvalues (sorted)');
cm = lines(run.NC_str);
figure;hold on;
% for k=1:run.NC_str
%     cells = find(run.labels_str == k);
%     scatter3(Y(cells,1),Y(cells,2),Y(cells,3),40,cm(k,:))
% end
for k=1:run.NC_str
    cells = find(run.labels_str == k);
    scatter(Y(cells,1),Y(cells,2),40,cm(k,:),'filled')
end
xlabel('1st component')
ylabel('2st component')



%% Four, but fourth one small
dim = sum(e > eps^(3/4))

% Poor reconstruction
maxerr2 = max(abs(pdist(X)-pdist(Y(:,1:2)))) 

% Good reconstruction
maxerr3 = max(abs(pdist(X)-pdist(Y(:,1:3)))) 

% Exact reconstruction
maxerr4 = max(abs(pdist(X)-pdist(Y)))

% D is now non-Euclidean
D = pdist(X,'cityblock');
[Y,e] = cmdscale(D);

% One is large negative
min(e)

% Poor reconstruction
maxerr = max(abs(pdist(X)-pdist(Y)))