function Points = CreateSobolNetwork(cellNo, meanSeperationDistance, SobolDimension)

initPointsNo = cellNo;
DistTolerance = 2;
p=sobolset(SobolDimension);
p = scramble(p,'MatousekAffineOwen');
Points = net(p,initPointsNo);

mDist = zeros((initPointsNo+1)^2,1);
for i=1:initPointsNo
    for j=1:initPointsNo
%         My mean is affected from zero values (diagonal)
            mDist((i-1)*(initPointsNo)+j) = mDist((i-1)*(initPointsNo)+j) + sqrt((Points(i,1)-Points(j,1))^2 + (Points(i,2)-Points(j,2))^2 + (Points(i,3)-Points(j,3))^2);
    end
end
% scrDist = mean(mDist);
% redoDist = (scrDist > meanSeperationDistance+DistTolerance) || (scrDist < meanSeperationDistance-DistTolerance);

% if redoDist
    scaleFactor = meanSeperationDistance / (1 / (cellNo^(1/3)));
%     DistFactor = scrDist / meanSeperationDistance;
%     Points = Points / DistFactor;
Points = Points * scaleFactor + (rand(1)*30);
% end


% fid = fopen('networkPositions.txt','w');
% for i=1:size(Points,1)
%     fprintf(fid, '%f\n',Points(i,1));
%     fprintf(fid, '%f\n',Points(i,2));
%     fprintf(fid, '%f\n',Points(i,3));
% end
% fclose(fid);


return;