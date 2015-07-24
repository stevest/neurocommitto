function [NCN] = m_commonNeighbors(adj)
    
NCN = zeros(size(adj));

adj = adj | adj' ;
% imagesc(adj);
% pause();

start=1;
for i=1:length(adj)
    for j=start:length(adj)   
        if (i~=j)
            NCN(i,j) = sum(all(adj(:,[i,j]),2));
        else
           NCN(i,j)=0;
        end
    end
%     start = start +1;
end


end