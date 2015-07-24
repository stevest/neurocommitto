function [CNI, CNO] = m_commonNeighbors(adj)
    
CNI = zeros(size(adj));
CNO = zeros(size(adj));

start=2;
for i=1:length(adj)
    for j=start:length(adj)
        CNI(i,j) = sum(all(adj(:,[i,j]),2));
        CNO(i,j) = sum(all(adj([i,j],:)));
    end
    start = start +1;
end

end