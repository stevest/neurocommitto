function Points = CreatePerinNetwork(step, jitter)
% Initializes neuronal somata as in example in Perin et al.
% First create a lattice of 12x12x13 with distances 'step' and afterwards
% add random special jitter with maximum value 'jitter'


ctr = 1;
for z=0:step:step*12
    for x=0:step:step*11
        for y=0:step:step*11
            Points(ctr,1:3) = ...
                [x+rand(1)*jitter,y+rand(1)*jitter,z+rand(1)*jitter] ;
            ctr = ctr+1;
        end
    end
end

return;