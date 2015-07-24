function equalizeClusters(cellsPerCluster,thrshld)

for d=1:length(cellsPerCluster)
    [~,tempIDX{d}] = sort((cellsPerCluster{d}),'descend');
end

ctr = 1;
while ~any(cellfun(@isempty,tempIDX))
    values = cellfun(@(x,idx) cellsPerCluster{idx}(x(1)),...
            tempIDX,num2cell(1:length(cellsPerCluster)));
    score(ctr) = sum(abs(diff(values)));
    [~,Idx2Remove(ctr)] = max(values);
    tempIDX{Idx2Remove(ctr)}(1) = [];
    ctr = ctr + 1;
end

%if we have score zero, great!
upTo = find(score==0);
if( upTo )
    for d=1:length(cellsPerCluster)
        [~,tempIDX{d}] = sort((cellsPerCluster{d}),'descend');
    end
    ctr = 1;
    while ctr<=upTo(1)
        tempIDX{Idx2Remove(ctr)}(1) = [];
        ctr = ctr + 1;
    end
else
end


end