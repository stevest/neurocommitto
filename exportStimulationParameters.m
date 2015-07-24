function exportStimulationParameters(AllCells,StimVect,outPath)
% Export stimulation parameters in .hoc file:

% Write clustering results to a .hoc file:
fid = fopen([outPath,'importStimulationParameters.hoc'],'w');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');
fprintf(fid,'// Override variables\n');

fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref PcellStimList\n');
fprintf(fid,'PcellStimList = new Vector(nPCcells)\n');

fprintf(fid,'\n\n// Import parameters:\n\n');
% network stimulation:
for i=1:AllCells(1)
    fprintf(fid,'PcellStimList.x[%d] = %d\n',i-1, StimVect(i));
end
fclose(fid);