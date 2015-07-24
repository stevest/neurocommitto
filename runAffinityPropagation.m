function [NC,labels,Sid] = runAffinityPropagation(mydata)
% Affinity Propagation clusteirng (see Frey & Dueck, Science, Feb. 2007)
% Note: Statistics Toolbox of Matlab needs to be installed

algorithm = 1;  % 1 --- adaptive AP, 0 --- original AP
nrun = 50000;   % max iteration times, default 50000
% nrun2 = 2000;   % max iteration times for original AP
nconv = 50;     % convergence condition, default 50
pstep = 0.01;   % decreasing step of preferences: pstep*pmedian, default 0.01
lam = 0.5;      % damping factor, default 0.5
cut = 3;        % after clustering, drop an cluster with number of samples < cut
%splot = 'plot'; % observing a clustering process when it is on
splot = 'noplot';


% initialization
type = 1;       % 1: Euclidean distances
simatrix = 1;   % 0: data as input; 1: similarity matrix as input
p = [];
Ms = [];
% derive similarity matrix
[rows,cols]=find(mydata);
for i=1:length(rows)
    M(i,1) = rows(i);
    M(i,2) = cols(i);
    M(i,3) = mydata(rows(i),cols(i));
end
data = [];
truelabels = ones(size(mydata,1),1);


disp(' '); disp(['==> Clustering is running, please wait ...']);
tic;
[labels,NCs,labelid,iend,Sp,Slam,NCfixs] = adapt_apcluster(M,type,...
    p,pstep,simatrix,'convits',nconv,'maxits',nrun,'dampfact',lam,splot);

[NC,Sil,Silmin] = solution_evaluation(data,M,labels,NCs,...
    NCfixs,simatrix,size(mydata,1),type,cut);
trun = toc;

fprintf('\n## Running time = %g seconds \n', trun);
fprintf('## Running iterations = %g \n', iend);

% finding an optimal clustering solution
solution_findK

truek = unique(truelabels);
truek = length(truek);
if truek > 1
    C = valid_external(labels(:,Sid), truelabels);
    fprintf('Fowlkes-Mallows validity index: %f\n', C(4));
end
if NCopt == truek
    fprintf('\n## Error rate of clustering solution might be inaccurate if large');
    fprintf('\n     (then use FM index instead) and it is for reference only:');
    valid_errorate(labels(:,Sid), truelabels);
end

end