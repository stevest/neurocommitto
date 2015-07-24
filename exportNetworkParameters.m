function exportNetworkParameters(AllCells,AllConnMat,AllWeightMat,outPath)
% Export parameter matrices in .hoc file:

% Write NMDA results to a .hoc file:
fid = fopen([outPath,'/importNetworkParameters.hoc'],'w');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');
fprintf(fid,'// Override variables\n');

fprintf(fid,sprintf('nPCcells=%d\n',AllCells(1)));
fprintf(fid,sprintf('nPVcells=%d\n',AllCells(2)));
fprintf(fid,sprintf('nCBcells=%d\n',AllCells(3)));
fprintf(fid,sprintf('nCRcells=%d\n',AllCells(4)));
fprintf(fid,sprintf('nAllCells=%d\n\n',AllCells(5)));

fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref connMatrix, weightsMatrix\n');
fprintf(fid,'connMatrix = new Matrix(nAllCells, nAllCells)\n');
fprintf(fid,'weightsMatrix = new Matrix(nAllCells, nAllCells)\n');

fprintf(fid,'\n\n// Import parameters: (long-long text following!)\n\n');
% network connectivity:
for i=1:length(AllConnMat)
    for j=1:length(AllConnMat)
        fprintf(fid,'connMatrix.x[%d][%d] = %d\n',i-1,j-1, AllConnMat(i,j));
    end
end
% Network synaptic weights
for i=1:length(AllConnMat)
    for j=1:length(AllConnMat)
        fprintf(fid,'weightsMatrix.x[%d][%d] = %f\n',i-1,j-1, AllWeightMat(i,j));
    end
end
fclose(fid);
end