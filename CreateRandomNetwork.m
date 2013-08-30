function Points = CreateRandomNetwork(cellNo, maxSeperationDistance, RandomDimension)

initPointsNo = cellNo;
DistTolerance = 2;
Points = rand(initPointsNo,RandomDimension);

% scatter3(Points(:,1), Points(:,2), Points(:,3));

mDist = zeros((initPointsNo+1)^2,1);
for i=1:initPointsNo
    for j=1:initPointsNo
%         My mean is affected from zero values (diagonal)
            mDist((i-1)*(initPointsNo)+j) = mDist((i-1)*(initPointsNo)+j) + sqrt((Points(i,1)-Points(j,1))^2 + (Points(i,2)-Points(j,2))^2 + (Points(i,3)-Points(j,3))^2);
    end
end
scrDist = max(mDist) ;
redoDist = (scrDist < maxSeperationDistance+DistTolerance) || (scrDist > maxSeperationDistance-DistTolerance);

if redoDist
    DistFactor = scrDist / maxSeperationDistance;
    Points = Points / DistFactor;
end


return;