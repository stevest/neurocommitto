function somata = CreateCubeNetworkPV(cubeSize, cellNo)
% Initializes neuronal PV somata in a cube of given dimentions
% arg_1 the value of cube side in ?m.
% arg_2 is the number of cells

initPointsNo = cellNo;
DistTolerance = 2;
somata = rand(initPointsNo,3);

somata = somata .* cubeSize;

return;