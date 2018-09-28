%% SCRIPT USED FOR MODEL GOODNESS-OF-FIT AND PARAMETER INFERENCE:
% THE log-LIKELIHOOD VALUE IS RETRIEVED & FITS ON CLONE SIZE DISTRIBUTIONS ARE PLOTTED

%% LOAD EXPERIMENTAL CLONE SIZE DATA AND CONVERT TO CLONE SIZE DISTRIBUTIONS:
load Lrig1_LineageTracing_dataset.mat
% rtime: vector containing experimental timepoints (expressed in weeks)
% rx_basal: cell array of format {nmice,timepoints} containing basal clone sizes per individual animal, individual time point
% rx_basal_all: cell array of format {1,timepoints} containing basal clone sizes per individual time point (animals pooled)

% Collect frequencies for each clone size (No. of basal cells):
% we may exclude clones with less than 2 basal cells, since in the experiments they may come from labelling postmitotic, already differentiating cells
[rfreq_all, rfreq_all_rel] = size2freq(rx_basal_all,rtime,1,rx_basal_all,2); % only clones with at least 2 basal cells are considered

% Present frequencies in bins (increasing in powers of 2):
[rfreq_bin_all, rfreq_bin_all_rel, rbin_label] = size2freqbinned(rfreq_all,rx_basal_all,rtime,1);

%% COMPUTATIONAL SIMULATION OF CLONE SIZES UNDER SPECIFIC PARAMETER CONDITIONS
% The single-progenitor (SP) model is simulated under specific parameter
% conditions and computational clone size distributions retrieved. They
% will constitute the PDF to be tested against experimental datasets

% Example with specific parameter values: (those from MLE on Lrig1 data)
% rtime: same as experimental time points
lambda = 2.9; %(/week)
r = 0.095;
gamma = 5.3857; %(/week)
indiv = 10000; % estimated calculation time: ~1min for 10k clones | ~10min for 100k clones
tlag = 0.5/7; % half a day
GamShape = 8;
[nx_basal,ntime] = MonteCarloSimulator_SP_BasalCloneDynamics(rtime,lambda,r,gamma,indiv,tlag,GamShape);

% Collect frequencies for each clone size (No. of basal cells):
% we may exclude clones with less than 2 basal cells if we did so with the experimental data sets
[nfreq, nfreq_rel] = size2freq(nx_basal,rtime,2,nx_basal,2);

% Present frequencies in bins (increasing in powers of 2):
[nfreq_bin, nfreq_bin_rel, nbin_label] = size2freqbinned(nfreq,nx_basal,rtime,2);

%% GOODNESS-OF-FIT CALCULATED AS log-LIKELIHOOD VALUE:
% A value of log-Likelihood is obtained for the specific set of parameter
% values simulated.
% This procedure can be repeated for multiple parameter values (if we were
% to loop on them, in a parameter sweep) to estimate the log-Likelihood
% landscape and infer the combination of values yielding the maximum
% likelihood estimate (MLE)

Lbasal_bin_t = [];
Lbasal_bin = [];
try
    %                                                       rfreq      myPDF     timepoints
    [Lbasal_bin_t(1,:), Lbasal_bin(1,1)] = logLike_calc(rfreq_bin_all,nfreq_bin_rel,rtime);
catch
    disp('Could not compute the LogLikelihood!!! Value assigned NaN')
    Lbasal_bin_t(1,:) = NaN(1,size(rtime,2)); Lbasal_bin(1,1) = NaN;
end

disp(sprintf('log-Likelihood value: %f',Lbasal_bin))

%% PLOT FITTING ON CLONE SIZE DISTRIBUTIONS OVER TIME:
figure()
mycol = [1 0 0; 0.93 0.69 0.13; 0.47 0.67 0.19; 0 0 1; 0.75 0 0.75; zeros(10,3)]; % color palette

% Experimental clone size distributions:
avg_rfreq_bin_all_rel = {}; sem_rfreq_bin_all_rel = {};
for aja = 1:length(rtime)
    % mean frequency of clones of certain size:
    avg_rfreq_bin_all_rel{1,aja} = rfreq_bin_all_rel{1,aja}(3:end,1) ./ sum(rfreq_bin_all_rel{1,aja}(3:end,1));
    % SEM of frequency of clones of certain size
    sem_rfreq_bin_all_rel{1,aja} = sqrt(avg_rfreq_bin_all_rel{1,aja} .* (1-avg_rfreq_bin_all_rel{1,aja})) ./ sqrt(sum(rfreq_bin_all{1,aja}(3:end,1)));
    % Plot as errorbars:
    hold on
    for lops = 1:size(avg_rfreq_bin_all_rel{1,aja},1)
        errorbar(rtime(aja),avg_rfreq_bin_all_rel{1,aja}(lops,1),sem_rfreq_bin_all_rel{1,aja}(lops,1),'o','Color',mycol(lops,:),'MarkerFaceColor',mycol(lops,:))
    end
end

% Computational clone size distributions:
avg_nfreq_bin_rel = [];
avg_nfreq_bin_rel = nfreq_bin_rel(3:end,:) ./ sum(nfreq_bin_rel(3:end,:),1);
% Plot frequency time courses:
for eje = 1:size(avg_nfreq_bin_rel,1)
    plot(ntime,avg_nfreq_bin_rel(eje,:),'-','Color',mycol(eje,:))
end

set(gca,'XScale','log'); box on; ylim([0 1]); %xlim([0.3 100]); 
xlabel('Time (weeks)'); ylabel('Clone frequency (%)')
%legend('2','3-4','5-8','9-16','17-32','33-64','>64')