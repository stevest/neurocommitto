function [within,between,WB] = calculateWithinBetween(conn,labels)

cellsPerCluster_str = {};
WB = {};
within = [];
between = [];
for k=1:size(conn,3)
    cellsPerCluster_str(k) = {histc(labels(:,k),1:length(unique(labels(:,k))))'};
    % calculate cluster inter/intra conn probs:
    for j = 1:size(cellsPerCluster_str{k},2)
        for l = 1:size(cellsPerCluster_str{k},2)
            jIdx = find(labels(:,k)==j);
            lIdx = find(labels(:,k)==l);
            % /!\ ATTENTION because: (A+B)/(Na+Nb) ~= (A/Na) + (B/Nb), we
            % do A/(Na+Nb) + B/(Na+Nb)
            % Also autapses do not count:
            if j ~= l
                WB{k}(j,l) = sum(sum(conn(jIdx,lIdx,k))) / (numel(conn(jIdx,lIdx,k))*2);
            else
                WB{k}(j,l) = sum(sum(conn(jIdx,lIdx,k))) / (numel(conn(jIdx,lIdx,k))-length(lIdx));
            end
        end
    end
    within(k) = nanmean(diag(WB{k}));
%     within(k) = nanmax(diag(WB{k}));
    %sum incomming/outgoing connections: 
    % The PLUS from the above ATTENTION:
    tmp = triu(WB{k} + WB{k}',1);
    tmp = tmp(logical(triu(ones(size(WB{k})),1)));
    if ~isempty(tmp)
    between(k) = nanmean(tmp);
%     between(k) = nanmax(tmp);
    else
        between(k) = NaN;
    end
end