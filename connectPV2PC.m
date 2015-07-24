function ConnMatPV2PC = connectPV2PC(dist,connProbsFunc)
% rng('shuffle');

ConnMatPV2PC = zeros(size(dist));
for i=1:size(dist,1)
    for j=1:size(dist,2)
        if(rand(1) <= connProbsFunc(dist(i,j)))
           ConnMatPV2PC(i,j) = 1;
        end
    end
end

return;
end
