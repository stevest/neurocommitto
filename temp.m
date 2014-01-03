nPC = 75;%216;
clstr = 1; % Cluster type
connBinsPC2PC = 20:30:500;
connProbsPC2PC = 0.25 .* exp(-0.006 * connBinsPC2PC);
connBinsPV2PC = 0:20:500;
connProbsPV2PC = linspace(1,0,length(connBinsPV2PC)) * 0.7;
connBinsCB2PC = 0:20:500;
connProbsCB2PC = normpdf(connBinsCB2PC,150,90) *10;
recipBinsPC2PC = 20:30:500;
recipProbsPC2PC = 0.12 .* exp(-0.006 * recipBinsPC2PC);
% Update probabilities from Perin et al. (Figure S4)
NconnBins = [0,1,2,3];
NincomingProbs = [0.1, 0.2, 0.25, 0.4] ;
NoutgoingProbs = [0.12, 0.25, 0.24, 0.2] ;
%  ---- Initialize Pyramidal cells ----
PCsomata = CreateRandomNetwork(nPC, 200, 3);
distPC2PC = generateDistanceMat(PCsomata', 0);
% Pyramidals connect to all
% glProbConn = 0.4;
PC2PC = zeros(nPC,nPC,4);
PC2PC(:,:,5) = create_graph_WS(distPC2PC,0.4,0.99); % 1.0=Random graph, 0.0=Watts-Strogatz graph
PC2PC(:,:,4) = create_graph_WS(distPC2PC,0.4,0.5); % 1.0=Random graph, 0.0=Watts-Strogatz graph
PC2PC(:,:,3) = create_graph_WS(distPC2PC,0.4,0.01); % 1.0=Random graph, 0.0=Watts-Strogatz graph
PC2PC(:,:,2) = create_graph_DD(distPC2PC,0.9,connBinsPC2PC, connProbsPC2PC);
PC2PC(:,:,1) = create_graph_CN(distPC2PC,0.9, connBinsPC2PC, connProbsPC2PC, ...
recipBinsPC2PC, recipProbsPC2PC,NconnBins,NincomingProbs,NoutgoingProbs);


strict(:,:,1) =  PC2PC(:,:,1) .*  PC2PC(:,:,1)';
loose(:,:,1) = PC2PC(:,:,1) |  PC2PC(:,:,1)';
strict(:,:,2) =  PC2PC(:,:,2) .*  PC2PC(:,:,2)';
loose(:,:,2) = PC2PC(:,:,2) |  PC2PC(:,:,2)';
strict(:,:,3) =  PC2PC(:,:,3) .*  PC2PC(:,:,3)';
loose(:,:,3) = PC2PC(:,:,3) |  PC2PC(:,:,3)';
strict(:,:,4) =  PC2PC(:,:,4) .*  PC2PC(:,:,4)';
loose(:,:,4) = PC2PC(:,:,4) |  PC2PC(:,:,4)';
strict(:,:,5) =  PC2PC(:,:,5) .*  PC2PC(:,:,5)';
loose(:,:,5) = PC2PC(:,:,5) |  PC2PC(:,:,5)';
loose = double(loose);

for i=2:5
    [StrictCliques1(i), blah] = Cliquer.FindAll(strict(:,:,1), i, i, false, 3);
%     [StrictCliques2(i), blah] = Cliquer.FindAll(strict(:,:,2), i, i, false, 3);
%     [StrictCliques3(i), blah] = Cliquer.FindAll(strict(:,:,3), i, i, false, 3);
%     [StrictCliques4(i), blah] = Cliquer.FindAll(strict(:,:,4), i, i, false, 3);
%     [StrictCliques5(i), blah] = Cliquer.FindAll(strict(:,:,5), i, i, false, 3);
%     
%     [LooseCliques1(i), blah] = Cliquer.FindAll(loose(:,:,1), i, i, false, 3);
%     [LooseCliques2(i), blah] = Cliquer.FindAll(loose(:,:,2), i, i, false, 3);
%     [LooseCliques3(i), blah] = Cliquer.FindAll(loose(:,:,3), i, i, false, 3);
%     [LooseCliques4(i), blah] = Cliquer.FindAll(loose(:,:,4), i, i, false, 3);
%     [LooseCliques5(i), blah] = Cliquer.FindAll(loose(:,:,5), i, i, false, 3);
end
