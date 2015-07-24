function ConnMatPV2PC = connectPV2PC(dist,connProbs)

ConnMatPV2PC = zeros(size(dist));

for i=1:size(dist,1)
    for j=1:size(dist,2)
        pp=connProbs(dist(i,j));
        if(rand(1) <= pp)
           ConnMatPV2PC(i,j) = 1;
        end
    end
end

return;
end
