function Points = CreateSliceNetwork(netSize, step, jitter)
% Initializes neuronal somata as in example in Perin et al.
% First create a lattice of 12x12x13 with distances 'step' and afterwards
% add random special jitter with maximum value 'jitter'

if netSize(3) == 0
    jitterZ = 0;
else
    jitterZ = rand(1);
end
rng('Shuffle');


ctr = 1;
for z=0:step:step*netSize(3)
    for x=0:step:step*netSize(2)
        for y=0:step:step*netSize(1)
            jitterX = rand(1);
            jitterY = rand(1);
            jitterVec = [jitterX, jitterY, jitterZ] ;
            jitterVec = (jitterVec/norm(jitterVec)) * jitter;
            Points(ctr,1:3) = ...
                [x+jitterVec(1),y+jitterVec(2),z+jitterVec(3)] ;
            ctr = ctr+1;
        end
    end
end

return;