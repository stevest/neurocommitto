function [Points, distPC2PC] = CreateRandomNetwork(cellNo, maxSeperationDistance)

initPointsNo = cellNo;
DistTolerance = 2;
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

for i=1:initPointsNo
    for j=1:initPointsNo
        distPC2PC(i,j) = sqrt((Points(i,1)-Points(j,1))^2 + (Points(i,2)-Points(j,2))^2 + (Points(i,3)-Points(j,3))^2); 
    end
end

return;