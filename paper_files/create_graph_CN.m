% Generate an instance of a distance dependent graph with a common neighbor
% rule. Consideres the distance & CN probabilities to be independent.
% Based (but modified) on Perin 2011/2013.
% Nassi, 9/03/15

function [GE,E, GCC,ACC, mean_NCN, pp]=create_graph_CN(distMat, Cprobs, Rprobs, Tprobs)

p_target=0.13;
N = length(distMat);
E = zeros(size(distMat));
NCN = zeros(size(distMat));

% Initialize network connecticity based on intersomatic distances.
% First create Pd matrix: for every pair the overall connection
% probability due to distance.

for i=1:N
    for j=1:N
        if (i~=j)
            Pd(i,j)=Tprobs(distMat(i,j));
        else
            Pd(i,j)=0;
        end
    end
end

% Update probabilities to target (0.13, Morishima Kawagushi)
Pd = (p_target/(sum(sum(Pd))/(numel(Pd)-N))) * Pd;

% Connect network with Pd (distance-dependent) probabilities.
start = 2;
for i=1:N
    for j=start:N
        xx=(log(0.22)-log(Pd(i,j)))/0.006;
        pnr = Cprobs(xx);
        pr  = Rprobs(xx);
        randomNum = rand(1);
        if(randomNum < pr)
            E(i,j) = 1;
            E(j,i) = 1;
        elseif ((randomNum >= pr) && (randomNum < pr+(pnr/2)))
            E(i,j) = 1;
        elseif ((randomNum >= pr+(pnr/2)) && (randomNum < (pr+pnr)))
            E(j,i) = 1;
        end
    end
    start = start+1;
end

% Initial Random Network
GE=E;

% for an arbitrary number of iterations (check that ACC/GCC converge).
loop_t=20;

pp(1)=sum(sum(E))/((N*N)-N);
[GCC(1),ACC(1), ~] = clust_coeff(E);
for t=1:loop_t
    % Use average & global clustering coefficient (Perin 2013, Figure S4 uses average).
    %     [GCC(t),ACC(t), ~] = clust_coeff(E);
    
    % Update connectivity matrix
    NCN= m_commonNeighbors(E);
    
    mean_NCN(t)=mean(NCN(~eye(N)));
    if (length(find(NCN==0))==0)
        t
    end
    
    CN_prob=((NCN)./max(NCN(~eye(N))));
    
    % Final probability matrix as the product of CN_p and Pd
    F_prob=CN_prob.*Pd;
    %     sum(sum(F_prob))/((N*N)-N)
    % Normalize the new probability so that all incoming
    % connections are made with prob as in the initial network
    F_prob=F_prob* p_target*((N*N)-N)/sum(sum(F_prob));
    
    E = zeros(size(distMat));
    start=2;
    for i=1:N
        for j=start:N
            % J is the source
            randomNum = rand(1);
            
            xx=(log(0.22)-log(F_prob(i,j)))/0.0052;
            pr= Rprobs(xx);
            pnr=Cprobs(xx);
            
            if randomNum<pr
                E(i,j)=1;
                E(j,i)=1;
            elseif ((randomNum >= pr) && (randomNum < pr+(pnr/2)))
                E(i,j) = 1;
            elseif ((randomNum >= pr+(pnr/2)) && (randomNum < (pr+pnr)))
                E(j,i) = 1;
            end
        end
        start=start+1;
    end
    
    if  t==loop_t
        %          figure; bar(histc(reshape(NCN, [numel(distMat),1]), 0:1:15));
        %          ylim([0 1500])
        [GCC(2),ACC(2), ~] = clust_coeff(E);
        pp(2)=sum(sum(E))/((N*N)-N);
        
    end
    
end