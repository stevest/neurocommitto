function distMat = distancePV2PC(PV,PC)
    distMat = zeros(size(PV,1),size(PC,1));
    for i=1:size(PV,1)
       distMat(i,1:size(PC,1)) = reshape(sqrt( sum((PC(:,1:3)-repmat(PV(i,1:3),[size(PC,1),1])).^2,2) ),1,[]);
    end
return;
end