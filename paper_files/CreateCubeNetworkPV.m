function [Points, distMat] = CreateCubeNetworkPV(cubeSize, cellNo,PCsomata )
% Initializes neuronal PV somata in a cube of given dimentions
% arg_1 the value of cube side in ?m.
% arg_2 is the number of cells

initPointsNo = cellNo;
DistTolerance = 2;
maxSeperationDistance=cubeSize;

Points = rand(initPointsNo,3);
% Points=sortrows(Points);
%Calculate distances
mDist = zeros((initPointsNo+1)^2,1);
for i=1:initPointsNo
    for j=1:initPointsNo
            %Calculate all distances
            mDist((i-1)*(initPointsNo)+j) = mDist((i-1)*(initPointsNo)+j) + sqrt((Points(i,1)-Points(j,1))^2 + (Points(i,2)-Points(j,2))^2 + (Points(i,3)-Points(j,3))^2);
    end
end
% Find max value
scrDist = max(mDist) ;
redoDist = (scrDist < maxSeperationDistance+DistTolerance) || (scrDist > maxSeperationDistance-DistTolerance);

%Reshape to max Sepearition Distance 
if redoDist
    DistFactor = scrDist / maxSeperationDistance;
    Points = Points / DistFactor;
end

% somata = somata .* cubeSize;

distMat = zeros(size(Points,1),size(PCsomata,1));

for i=1:size(Points,1)
   distMat(i,1:size(PCsomata,1)) = reshape(sqrt(sum((PCsomata(:,1:3)-repmat(Points(i,1:3),[size(PCsomata,1),1])).^2,2) ),1,[]);
end


return;