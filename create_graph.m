% Generate an instance of a random graph

function [E]=create_graph(distMat, param)
N = length(distMat);
E = zeros(size(distMat));
for i=1:N % for each node:
   [~,IX] = sort(distMat(i,:),'ascend') ;
   for el=1:N/2 % connect about half the nodes together
       E(i,IX(el)) = 1;
   end
end

% start = 2;
for i=1:N
   for j=1:N
       % with random propability rellocate some connections
       if( (rand(1)<param) && E(i,j) )
           E(i,j) = 0;
           ri = round(rand(1,2)*(N-1))+1 ;
           while(E(ri(1),ri(2)))
               ri = round(rand(1,2)*(N-1))+1 ;
           end
           E(ri(1),ri(2)) = 1;
       end
   end
%    start = start+1;   
end
   
end 