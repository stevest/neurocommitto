function ConnMatCB2PC = connectCB2PC(dist,connBins,connProbs)
% rng('shuffle');

ConnMatCB2PC = zeros(size(dist));
for i=1:size(dist,1)
    for j=1:size(dist,2)
        [~,bin] = histc(dist(i,j),connBins);
        if(rand(1) <= connProbs(bin))
           ConnMatCB2PC(i,j) = 1;
        end
    end
end

return;
end
