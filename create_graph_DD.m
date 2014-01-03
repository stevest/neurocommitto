% Generate an instance of a distance dependent graph

function [E]=create_graph_DD(distMat,prob, bins, probs)
N = length(distMat);
E = zeros(size(distMat));

for i=1:N % for each node:
    for j=1:N
        tempIDX = find(histc(distMat(i,j),bins)) ;
        if(~isempty(tempIDX)) && (prob>rand(1))
            if(probs(tempIDX) >rand(1))
                E(i,j) = 1;
            end
        end
    end
end

end