function distMat = distancePV2PC(PV,PC)
    distMat = zeros(length(PV),length(PC));
    for i=1:length(PV)
       distMat(i,1:length(PC)) = sqrt( sum((PC(:,1:3)-repmat(PV(i,1:3),[length(PC),1])).^2,2) );
    end
return;
end