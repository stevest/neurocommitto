function exportNetworkInhibition(InhibInputGABAaSoma,outPath)
% Export GABAa inhibition (Background) parameters in .hoc file:
fid = fopen([outPath,'/importNetworkInhibition.hoc'],'w');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');

fprintf(fid,'// Object decleration:\n');
fprintf(fid,'BG_somaSyn = %d\n', size(InhibInputGABAaSoma,2));
fprintf(fid,'objref Inhib_GABAa_Soma[%d][%d]\n',...
    length(InhibInputGABAaSoma),size(InhibInputGABAaSoma,2));

fprintf(fid,'\n\n// Import parameters: \n\n');
% Only for Pyramidals:
for c=1:size(InhibInputGABAaSoma,1)
    for i=1:size(InhibInputGABAaSoma,2)
        fprintf(fid,'Inhib_GABAa_Soma[%d][%d] = new Vector(%d)\n',c-1,i-1,length(InhibInputGABAaSoma{c,i}));
        for j=1:length(InhibInputGABAaSoma{c,i})
            fprintf(fid,'Inhib_GABAa_Soma[%d][%d].x[%d] = %d\n',c-1,i-1,j-1, abs(InhibInputGABAaSoma{c,i}(j)));
        end
    end
end

fclose(fid);
end
