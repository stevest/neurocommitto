function rndCombs = combntnsRND(n,k,maxCombs)
% Produces all OR a number of possible K combinations from a set N.
% If all possible combinations are few (can be handled by hardware),
% produces all the possible combinations as in :
% http://en.wikipedia.org/wiki/Combination
% using the Matlab's build in 'combntns' function.
% Else produces a random numder of combinations in an efford to cover the
% combination space as even as possible (but NOT all possible combinations
% are concidered..)
%
% stamatiad.st@gmail.com

feasible = factorial(n) / ( factorial(k) * factorial(n-k) );

if(~isnan(feasible)) && (feasible ~= inf)
    rndCombs = combntns(1:n,k);
    return;
else
    rndCombs=[];
    i=1;
    while i < maxCombs
        randK = ceil(rand(1,k) * n) ;
        if ( ~any(ismember(randK, rndCombs,'rows')) )
            rndCombs = [rndCombs; randK] ;
            i = i+1;
        end
    end
    return;
end
return;

end