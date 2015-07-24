% clc;clear all;close all;
% create dummy data:
% ST_1 = find(rand(1,18000)>0.5);
% ST_2 = find(rand(1,18000)>0.5);
nIncoming=[250,250]
randShift = 50%71
nSpikes = 20
tstop = 5000
spikesDend={};
spikesApic={};
for c=1:75
    i=1;
    while i <= nIncoming(1)
        spikesDend(c,i)={(rand(1)*500) + linspace(0,tstop,nSpikes)}; % random phase
        spikesDend{c,i} = round(spikesDend{c,i} + ((rand(1,nSpikes)-0.5)*randShift));
        spikesDend{c,i} = spikesDend{c,i}(spikesDend{c,i}>0);
        spikesDend{c,i} = spikesDend{c,i}(spikesDend{c,i}<tstop);
        spikesDend{c,i} = sort((spikesDend{c,i}),'ascend');
        if(length(unique(spikesDend{c,i})) ~= length(spikesDend{c,i}) )
            continue;
        end
        i = i +1;
    end
    i=1;
    while i <= nIncoming(2)
        spikesApic(c,i)={(rand(1)*500) + linspace(0,tstop,nSpikes)};
        spikesApic{c,i} = round(spikesApic{c,i} + ((rand(1,nSpikes)-0.5)*randShift));
        spikesApic{c,i} = spikesApic{c,i}(spikesApic{c,i}>0);
        spikesApic{c,i} = spikesApic{c,i}(spikesApic{c,i}<tstop);
        spikesApic{c,i} = sort((spikesApic{c,i}),'ascend');
        if(length(unique(spikesApic{c,i})) ~= length(spikesApic{c,i}) )
            continue;
        end
        i = i +1;
    end
end

nIncomingS=[60]
nSpikesS = 50 % ~6 Hz
tstop = 5000
randShiftS = tstop / nSpikesS * 1000
spikesSoma={};
for c=1:75
    i=1;
    while i <= nIncomingS(1)
        spikesSoma(c,i)={linspace(0,tstop,nSpikesS)};
        spikesSoma{c,i} = round(spikesSoma{c,i} + ((rand(1,nSpikesS)-0.5)*randShiftS));
        spikesSoma{c,i} = spikesSoma{c,i}(spikesSoma{c,i}>0);
        spikesSoma{c,i} = spikesSoma{c,i}(spikesSoma{c,i}<tstop);
        spikesSoma{c,i} = unique(sort((spikesSoma{c,i}),'ascend'));
%         if(length(unique(spikesSoma{c,i})) ~= length(spikesSoma{c,i}) )
%             continue;
%         end
        i = i +1;
    end
end


%%
% create dummy data:
% ST_1 = find(rand(1,18000)>0.5);
% ST_2 = find(rand(1,18000)>0.5);
nIncoming=[30,30]
 % Background activity must result in ~1.2Hz in the post-synaptic cell as in
 % Carl C.H.Petersen Neuron, 2013
tstop = 5000
spikesDend={};
spikesApic={};
randShift = 180;
randShiftS = 180;
stimPoisson = find(poissrnd(0.08,tstop,1))  ;
poil = length(stimPoisson) ;
for c=1:nPC
    for i = 1:nIncoming(1)
    spikesDend(c,i) = { sort(round((rand(poil,1)-0.5)*randShift) + stimPoisson,'ascend')} ; %{ find(poissrnd(0.004,tstop,1)) } ; % ~12Hz
    spikesDend{c,i}(spikesDend{c,i}<0) = [];
    spikesDend{c,i}(spikesDend{c,i}>tstop) = [];
    end

    for i = 1:nIncoming(2)
    spikesApic(c,i) = { sort(round((rand(poil,1)-0.5)*randShift) + stimPoisson, 'ascend')}  ; % ~12Hz
    spikesApic{c,i}(spikesApic{c,i}<0) = [];
    spikesApic{c,i}(spikesApic{c,i}>tstop) = [];
    end
end
% symfwna me kapoion Calr Petersen prepei to spiking tou inhibition (~11Hz)na
%   spiking (freq) of excitation (~0.35Hz)
nIncomingS=[4]
spikesSoma={};
for c=1:nPV
    for i = 1:nIncomingS(1)
    spikesSoma(c,i) = { sort(round((rand(poil,1)-0.5)*randShiftS) + stimPoisson,'ascend')} ; % ~12Hz
    spikesSoma{c,i}(spikesSoma{c,i}<0) = [];
    spikesSoma{c,i}(spikesSoma{c,i}>tstop) = [];
    end
end

%%
mypath = '/home/stefanos/Documents/prefrontal-micro/experiment/network/';

exportBackgroundStimParams(spikesDend,spikesApic,spikesSoma,mypath);


%%
spikesGeneral = spikesDend ;
% Init SPIKY:
para.tmin = 0;
para.tmax =  5000;
para.dts =  10;
para.dtm = 10;
para.select_measures=[0 1     0 0     0 0     0 0];            % Select order of measures


% feed SPIKY

results = SPIKY_no_plot_f_distances_MEX(spikesGeneral(1,:),para);

% ################################################### Example plotting (just as a demonstration)
plotting=7;
num_plots=(mod(plotting,2)>0)+(mod(plotting,4)>1)+(mod(plotting,8)>3);

measures=find(para.select_measures);
for mc=1:length(measures)
    measure=measures(mc);
    figure(mc); clf
    set(gcf,'Name',results.measures{mc})
    set(gcf,'Units','normalized','Position',[0.0525 0.0342 0.8854 0.8867])
    subplotc=0;
    
    if mod(plotting,2)>0
        subplotc=subplotc+1;
        subplot(num_plots,1,subplotc)                                                                     % Spikes
        for trc=1:length(spikesGeneral)
            for spc=1:length(spikesGeneral{trc})
                line(spikesGeneral{trc}(spc)*ones(1,2),[trc-1 trc])
            end
        end
        xlim([para.tmin para.tmax])
        ylim([0 length(spikesGeneral)])
        title ('Spike trains','FontWeight','bold','FontSize',14)
    end

    if mod(plotting,4)>1
        subplotc=subplotc+1;
        subplot(num_plots,1,subplotc)                                                                     % Dissimilarity profile
        if ismember(measure,[1 3 4])                 % PICO profiles have first to be transformed
            [overall_dissimilarity,plot_x_values,plot_y_values] = SPIKY_f_pico(results.isi,results.dissimilarity_profiles{mc},para.dts,para.tmin);
            plot(plot_x_values,plot_y_values)
        elseif ismember(measure,[2 7 8])             % piecewise linear profiles can be plotted right away                                                                % sampled profiles can be plotted right away
            plot(results.pili_time,results.dissimilarity_profiles{mc})
        else                                                                % sampled profiles can be plotted right away
            plot(results.time,results.dissimilarity_profiles{mc})
        end
        xlim([para.tmin para.tmax])
        title ([results.measures{mc},'   ---   Dissimilarity profile'],'FontWeight','bold','FontSize',14)
    end

    if mod(plotting,8)>3
        subplotc=subplotc+1;
        subplot(num_plots,1,subplotc)                                                                      % Dissimilarity matrix
        if length(size(results.distance_matrices))==2
            mat=results.distance_matrices;
        else
            mat=shiftdim(results.distance_matrices(mc,:,:),1);
        end
        imagesc(mat)
        axis square
        title ([results.measures{mc},'   ---   Dissimilarity matrix'],'FontWeight','bold','FontSize',14)
    end
end
%%
signal_r = zeros(9000,1) ;
signal_r(round(rand(90,1) * 9000)) = 1;

signal_p = poissrnd(0.01,9000,1) ;


m = length(signal_r);          % Window length
n = pow2(nextpow2(m))
y=fft(signal_r,n);
f = (0:n-1);     % Frequency range
power = y.*conj(y)/n;   % Power of the DFT

plot(f,power)
xlabel('Frequency (Hz)')
ylabel('Power')
title('{\bf Periodogram}')
%%
BIG = zeros(1,5000);
for i= 1:size(spikesDend,1) * size(spikesDend,2)
    BIG(spikesDend{i}) = BIG(spikesDend{i}) + 1;
end