function exportNetworkStimulation(StimInputDend,StimInputApic,outPath)
% Export stimulation parameters in .hoc file:
fid = fopen([outPath,'/importNetworkStimulation.hoc'],'w');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');

fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref Stim_Dend[%d][%d], Stim_Apic[%d][%d]\n',...
    length(StimInputDend),length(StimInputApic),size(StimInputDend,2),size(StimInputApic,2));

fprintf(fid,'\n\n// Import parameters: \n\n');
% Only for Pyramidals:
for c=1:size(StimInputDend,1)
    for i=1:size(StimInputDend,2)
        fprintf(fid,'Stim_Dend[%d][%d] = new Vector(%d)\n',c-1,i-1,length(StimInputDend{c,i}));
        for j=1:length(StimInputDend{c,i})
            fprintf(fid,'Stim_Dend[%d][%d].x[%d] = %d\n',c-1,i-1,j-1, abs(StimInputDend{c,i}(j)));
        end
    end
end


for c=1:size(StimInputApic,1)
    for i=1:size(StimInputApic,2)
        fprintf(fid,'Stim_Apic[%d][%d] = new Vector(%d)\n',c-1,i-1,length(StimInputApic{c,i}));
        for j=1:length(StimInputApic{c,i})
            fprintf(fid,'Stim_Apic[%d][%d].x[%d] = %d\n',c-1,i-1,j-1, abs(StimInputApic{c,i}(j)));
        end
    end
end

fclose(fid);
end
