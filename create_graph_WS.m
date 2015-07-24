% Generate an instance of a random graph

function [E]=create_graph_WS(distMat, prob,flag)
N = length(distMat);
E = zeros(size(distMat));
% for i=1:N % for each node:
%     for el=1:N
%         if(rand(1) < prob)
%             E(i,el) = 1;
%         end
%     end
% end

% if strcmp(flag,'joint')
%     % As in Song et al, 2005
%     RND = rand( (numel(distMat)-size(distMat,1))/2,1 );
%     AB_IDX = findTriIdx(RND < prob, size(distMat,1),0) ;
%     E( AB_IDX ) = 1;
%     BA_IDX = findTriIdx((RND >= prob)&(RND<prob*2), size(distMat,1),1) ;
%     E( BA_IDX ) = 1;
%     RE_IDX = findTriIdx((RND >= prob*2)&(RND<(prob*2)+(prob^2)), size(distMat,1),0) ;
%     E( RE_IDX ) = 1;
%     RE_IDX = findTriIdx((RND >= prob*2)&(RND<(prob*2)+(prob^2)), size(distMat,1),1) ;
%     E( RE_IDX ) = 1;
% end




% if strcmp(flag,'ind')
%     % As in Song et al, 2005
%     RND = rand( numel(distMat) - size(distMat,1),1 );
%     AB_IDX = find(RND < (prob*(1-prob)) ) ; % 1436
%     E( AB_IDX ) = 1;
%     BA_IDX = find((RND >= (prob*(1-prob))) & (RND<(prob*(1-prob)*2)) ) ; % 1343
%     E( BA_IDX ) = 1;
% %     RE_IDX = find((RND >= (prob*(1-prob)*2))&(RND<((prob*(1-prob)*2)+(prob^2)))) ;
% %     E( RE_IDX ) = 1;
%     RE_IDX = find((RND >= (prob*(1-prob)*2))&(RND<(prob^2))) ;
%     E( RE_IDX ) = 1;
% end

% if strcmp(flag, 'ind')
%     % connect with independent probabilities:
%     E = (rand(size(distMat)) < prob);
% end

if strcmp(flag,'ind')
    % As in Song et al, 2005
    RND = rand( size(distMat) );
    E( RND < prob ) = 1;
    E(find(eye(length(distMat)))) = 0;
end

% start = 2;
% for i=1:N
%     for j=start:N
%         % with random propability rellocate some connections
%         if( (rand(1)<param) && (E(i,j)))
%             ri = round(rand(1,2)*(N-1))+1 ;
%             E(ri(1),ri(2)) = E(i,j);
%         end
%     end
%     start = start+1;
% end

end

% % Generate an instance of a random graph
% function [E]=create_graph_WS(distMat, prob, param)
% N = length(distMat);
% E = zeros(size(distMat));
% for i=1:N % for each node:
%     [~,IX] = sort(distMat(i,:),'ascend') ;
%     for el=1:N
%         if(rand(1) < prob)
%             E(i,IX(el)) = 1;
%         end
%     end
% end
%
% % start = 2;
% for i=1:N
%     for j=1:N
%         % with random propability rellocate some connections
%         if( (rand(1)<param) && (E(i,j)))
%             E(i,j) = 0;
%             ri = round(rand(1,2)*(N-1))+1 ;
%             while(E(ri(1),ri(2)))
%                 ri = round(rand(1,2)*(N-1))+1 ;
%             end
%             E(ri(1),ri(2)) = 1;
%         end
%     end
% %     start = start+1;
% end
%
% end

function outpt = findTriIdx(inpt, sizeof, flip)
outpt = [];
idx = find(inpt);
if numel(idx) == 0
    return;
end
outpt = zeros(length(idx),1);

ctr = 1;
arrctr = 1;
start = 2;
for i=1:sizeof
    for j=start:sizeof
        if(idx(ctr) == arrctr)
            if flip
                outpt(ctr) = sub2ind([sizeof,sizeof], j,i) ;
            else
                outpt(ctr) = sub2ind([sizeof,sizeof], i,j) ;
            end
            if ctr >=length(idx)
                return;
            end
            ctr = ctr+1;
        end
        arrctr = arrctr+1;
    end
    start = start+1;
end
end