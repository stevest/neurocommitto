function ConnMatPV2PC = connectPV2PC(dist,connBins,connProbs)
% rng('shuffle');

ConnMatPV2PC = zeros(size(dist));
for i=1:size(dist,1)
    for j=1:size(dist,2)
        [~,bin] = histc(143,connBins);
        if(rand(1) <= connProbs(bin))
           ConnMatPV2PC(i,j) = 1;
        end
    end
end

return;
end
