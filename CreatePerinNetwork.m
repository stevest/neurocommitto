function Points = CreatePerinNetwork(step, jitter)

ctr = 1;
for x=0:step:step*6%12
    for y=0:step:step*6%12
        for z=0:step:step*6%11
            Points(ctr,1:3) = ...
                [x+rand(1)*jitter,y+rand(1)*jitter,z+rand(1)*jitter] ;
            ctr = ctr+1;
        end
    end
end

return;