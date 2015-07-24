% Generate an instance of a distance dependent graph

function [E]=create_graph_DD(distMat, bins, probs)
% N = length(distMat);
E = zeros(size(distMat));


dist_probs = arrayfun(@(x) find(histc(x,bins)),distMat,'uniformoutput',false);
dist_probs = cellfun(@(x) probs(x), dist_probs,'uniformoutput',false) ;
dist_probs(cellfun(@isempty, dist_probs)) = {0};
dist_probs = cell2mat(dist_probs);
dist_probs = dist_probs(logical(tril(ones(size(dist_probs)),-1)));


% As in Song et al, 2005
RND = rand( length(dist_probs),1 );
AB_IDX = findTriIdx(RND < dist_probs, size(distMat,1),0) ;
E( AB_IDX ) = 1;
BA_IDX = findTriIdx((RND >= dist_probs)&(RND<dist_probs*2), size(distMat,1),1) ;
E( BA_IDX ) = 1;
RE_IDX = findTriIdx((RND >= dist_probs*2)&(RND<(dist_probs*2)+(dist_probs.^2)), size(distMat,1),0) ;
E( RE_IDX ) = 1;
RE_IDX = findTriIdx((RND >= dist_probs*2)&(RND<(dist_probs*2)+(dist_probs.^2)), size(distMat,1),1) ;
E( RE_IDX ) = 1;




% for i=1:N % for each node:
%     for j=1:N
%         tempIDX = find(histc(distMat(i,j),bins)) ;
%         if(~isempty(tempIDX)) && (prob>rand(1))
%             if(probs(tempIDX) >rand(1))
%                 E(i,j) = 1;
%             end
%         end
%     end
% end

end


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