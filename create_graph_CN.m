% Generate an instance of a distance dependent graph

function [E,CC]=create_graph_CN(distMat, Cbins, Cprobs, Rbins, Rprobs,Nbins, NIprobs, NOprobs)
N = length(distMat);
E = zeros(size(distMat));
CNI = zeros(size(distMat));
CNO = zeros(size(distMat));

% initialize based on intersomatic distance:
start = 2;
for i=1:N % for each node:
    for j=start:N
        % DOES NOT COUNT VALUES BIGGER THAN LAST VALUE OF BINS!!!
        tmpR = find(histc(distMat(i,j),Rbins)) ;
        tmpNR = find(histc(distMat(i,j),Cbins)) ;
        if( ~isempty(tmpR) && ~isempty(tmpNR)) 
            pnr = Cprobs(tmpNR);
            pr  = Rprobs(tmpR);
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
    end
    start = start+1;
end


% for t=1:10
% % sprintf('Overall ClustCoeff = %f\n',clust_coeff(E))
% CC(t) = clust_coeff(E);
%     % Update it
%     [CNI , CNO]= m_commonNeighbors(E);
%     E = zeros(size(distMat));
%     CNI((CNI>4)) = 4; % no more than 4 common neighbors
%     CNO((CNO>4)) = 4; % no more than 4 common neighbors
%     
%     start = 2;
%     for i=1:N % for each node:
%         for j=start:N
%             % DOES NOT COUNT VALUES BIGGER THAN LAST VALUE OF BINS!!!
%             tmpR = find(histc(distMat(i,j),Rbins)) ;
%             tmpNR = find(histc(distMat(i,j),Cbins)) ;
%             tmpCNIP = find(histc(CNI(i,j),Nbins)) ;
%             tmpCNOP = find(histc(CNO(i,j),Nbins)) ;
%             if( ~isempty(tmpR) && ~isempty(tmpCNIP) && ~isempty(tmpCNOP) )
%                 
%                 p = max([NIprobs(tmpCNIP), NOprobs(tmpCNOP)]);
%                 pnr = p*2;
% %                 pnr = Cprobs(tmpNR) * p ;
% %                 pr  = Rprobs(tmpR) *p ;%*0.8;
% %                 pnr = p;
%                 pr = pnr^2;
%                 if ((pnr/2)+pr) > ( (Cprobs(tmpNR)/2) + Rprobs(tmpR) )
% %                     sprintf('Overall Conn Prob = %f\n',Cprobs(tmpNR))
%                     f = ( (Cprobs(tmpNR)/2) + Rprobs(tmpR) ) / ((pnr/2)+pr) ;
%                     pnr = pnr*f;
%                     pr = pr*f;
%                 end
% 
% 
% %                 pr  = Rprobs(tmpR);%*0.8;
% %                 p = (NIprobs(tmpCNIP) + NOprobs(tmpCNOP));% * 0.08; %sss
% %                 pnr = (p - pr)*2; % Why not divided?
% % %                   p = (NIprobs(tmpCNIP) + NOprobs(tmpCNOP));
% % %                   pr = Rprobs(tmpR)*(p)*2;
% % %                   pnr = (p);
% % % %                 pnr = (p - pr);% my thought
% 
%                 randomNum = rand(1);
%                 if(randomNum < pr)
%                     E(i,j) = 1;
%                     E(j,i) = 1;
%                 elseif ((randomNum >= pr) && (randomNum < pr+pnr/2))
%                     E(i,j) = 1;
%                 elseif ((randomNum >= pr+pnr/2) && (randomNum < pr+pnr))
%                     E(j,i) = 1;
%                 end
%             else disp('FUCK!')
%             end
%         end
%         start = start+1;
%     end
%     
% end

end