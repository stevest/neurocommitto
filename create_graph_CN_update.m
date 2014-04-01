% Generate an instance of a distance dependent graph

function [E]=create_graph_CN_update(ConnMat,distMat, Cbins, Cprobs, Rbins, Rprobs,Nbins, NIprobs, NOprobs)
N = length(distMat);
E = zeros(size(distMat));
CNI = zeros(size(distMat));
CNO = zeros(size(distMat));

for t=1:10
    % Update it
    [CNI , CNO]= commonNeighbors(ConnMat);
    CNI(find(CNI>4)) = 4; % no more than 4 common neighbors
    CNO(find(CNO>4)) = 4; % no more than 4 common neighbors
    
    start = 2;
    for i=1:N % for each node:
        for j=start:N
            % DOES NOT COUNT VALUES BIGGER THAN LAST VALUE OF BINS!!!
            tmpR = find(histc(distMat(i,j),Rbins)) ;
            tmpCNIP = find(histc(CNI(i,j),Nbins)) ;
            tmpCNOP = find(histc(CNO(i,j),Nbins)) ;
            if( ~isempty(tmpR) && ~isempty(tmpCNIP) && ~isempty(tmpCNOP) )
                pr  = Rprobs(tmpR);
                p = NIprobs(tmpCNIP) + NOprobs(tmpCNOP);
                pnr = (p - pr)*2;
                randomNum = rand(1);
                if(randomNum < pr)
                    E(i,j) = 1;
                    E(j,i) = 1;
                elseif ((randomNum >= pr) && (randomNum < pr+pnr/2))
                    E(i,j) = 1;
                elseif ((randomNum >= pr+pnr/2) && (randomNum < pr+pnr))
                    E(j,i) = 1;
                end
            end
        end
        start = start+1;
    end
    
end

end