%% LOAD EXPERIMENTAL DATA OF H2BGFP INTENSITIES IN OESOPHAGEAL BASAL CELL NUCLEI AT DIFFERENT TIMES:
load H2BGFP_Int_OE_dataset.mat

%% LOG-2 TRANSFORM EXPERIMENTAL DATA, NORMALIZE IT, AND CALCULATE AVERAGE DIVISION RATE (for convenience):
% Initial definition of parameters:
log2_H2BGFP_Int = {};
log2n_H2BGFP_Int = {};
log2n_H2BGFP_Int_all = {};
avg_log2_H2BGFP_Int = [];
rtime_reshaped = [];

% Data renormalization:
for aja = 1:length(rtime)
    temporal_avg_log2_H2BGFP_Int = [];
    for eje = 1:nmice(aja)
        log2_H2BGFP_Int{eje,aja} = log2(H2BGFP_Int{eje,aja});
        temporal_avg_log2_H2BGFP_Int = [temporal_avg_log2_H2BGFP_Int mean(log2_H2BGFP_Int{eje,aja})];
        rtime_reshaped = [rtime_reshaped rtime(aja)];
    end
    temporal_renormalized_data = [];
    for eje = 1:nmice(aja)
        log2n_H2BGFP_Int{eje,aja} = log2_H2BGFP_Int{eje,aja} - mean(log2_H2BGFP_Int{eje,aja}) + mean(temporal_avg_log2_H2BGFP_Int);
        temporal_renormalized_data = [temporal_renormalized_data; log2n_H2BGFP_Int{eje,aja}];
    end
    log2n_H2BGFP_Int_all{1,aja} = temporal_renormalized_data;
    
    avg_log2_H2BGFP_Int = [avg_log2_H2BGFP_Int temporal_avg_log2_H2BGFP_Int];
end

% Calculate average division rate (from the linear slope of the log2(H2BGFP) decay)
myfun = fittype('p1*x + p2','dependent',{'y'},'independent',{'x'},'coefficients',{'p1', 'p2'});
[LineFit_all,LineFit_all_stat] = fit(rtime_reshaped',avg_log2_H2BGFP_Int',myfun,'Startpoint',[0 0],'Lower',[-Inf -Inf],'Upper',[Inf Inf]);
lambda_avg = abs(LineFit_all.p1);
disp(sprintf('Avg. Division rate: lambda = %.2f/week',lambda_avg));

%% Run ABC rejection method to deduce cell-cycle period distributions compatible with the actual pattern of H2BGFP dilution:
% Initialization parameters:
N = 10; % No. of acceptable posterior estimates
TolThres = 0.15; % Tolerance threshold for parameter acceptance/rejection
lambda = lambda_avg; % Range of possible values for the division rate
tlag_range = [0:0.25:2]./7; % Range of possible values for the refractory period (minimum cell-cycle period)
GamShape_range = 2.^[0:6]; % Range of possible values for the 'Shape' parameter of the gamma-distributed cell-cycle period

% Run algorithm:
[OK_lambda,OK_tlag,OK_GamShape] = ABCrejection_tccDist(rtime,log2n_H2BGFP_Int_all,N,TolThres,lambda,tlag_range,GamShape_range);

%% PLOT HEATMAP OF MOST-LIKELY tcc PROPERTIES (RefractPeriod & GamShape accepted parameter combinations):
figure()
OKmatrix = zeros((size(tlag_range,2)-1),size(GamShape_range,2));
for aja = 1:N
    ycoord = find(GamShape_range == OK_GamShape(aja));
    mytlag_pos = find(OK_tlag(aja) < tlag_range);
    OKmatrix(mytlag_pos(1)-1,ycoord) = OKmatrix(mytlag_pos(1)-1,ycoord) + 1;
end
imagesc(OKmatrix)
mycolormap = colormap('gray'); colormap(1-mycolormap)
set(gca,'XTick',1:7); set(gca,'YTick',1:8);
set(gca,'XTickLabel',num2cell(GamShape_range)); set(gca,'YTickLabel',{'0.0 - 0.25','0.25 - 0.5','0.5 - 0.75','0.75 - 1.0','1.0 - 1.25','1.25 - 1.5','1.5 - 1.75','1.75 - 2.0'});
ylabel('Refractory period (days)'); xlabel('Cell-cycle \Gamma-distribution shape, k');

%% PLOT ACCEPTED (DECONVOLUTED) CELL-CYCLE PERIOD DISTRIBUTIONS:
figure()
hold on
% first accepted cell-cycle period distributions fitting experimental observations:
for aja = 1:min([50 N])
    plot(([0:0.01:10]+OK_tlag(aja)).*7,gampdf([0:0.01:10],OK_GamShape(aja),(1/OK_lambda(aja)-OK_tlag(aja))./OK_GamShape(aja))./max(gampdf([0:0.01:10],OK_GamShape(aja),(1/OK_lambda(aja)-OK_tlag(aja))./OK_GamShape(aja))),'Color',[0.8 0.8 0.8]); xlim([0 14]);
end
% average cell-cycle period time:
line([1/mean(OK_lambda) 1/mean(OK_lambda)].*7,[0 1],'Color','k');
% example of an adequate cell-cycle period distribution:
mytlag = 0.5/7; myGam = 8;
plot(([0:0.01:10]+mytlag).*7,gampdf([0:0.01:10],myGam,(1/lambda_avg-mytlag)./myGam)./max(gampdf([0:0.01:10],myGam,(1/lambda_avg-mytlag)./myGam)),'Color','g'); xlim([0 14]);
hold off
ylabel('Frequency'); xlabel('Cell-cycle time (days)')
xlim([0 10]); set(gca,'XTick',[0:2:10]); set(gca,'YTick',[0:0.2:1]);

%% PLOT EXAMPLE OF FITS ON EXPERIMENTAL H2BGFP DILUTION HISTOGRAMS:
% Initialization parameters:
M = 10000; % No. of simulated individual basal cells per trial
tlag_all = [0, 0.5]./7;%2/7
GamShape_all = 2.^[0,3];%64
mystyle = {'--','-'};

figure()
for sj = 1:length(tlag_all)
    
    lambda = lambda_avg;
    tlag = tlag_all(sj);
    GamShape = GamShape_all(sj);

    % Simulation of H2BGFP dilution time course:
    Idist = MonteCarloSimulator_BasalCell_H2BGFPdil(rtime,log2n_H2BGFP_Int_all{1,1},lambda,M,tlag,GamShape);
        
    % Simulated H2BGFP intensity values are translated into log2 values and these normalized to avg at time0 (as in experimental data)
    clog2_Idist = log2(Idist)-median(log2(Idist(1,:)),2); % simulated
    rlog2_Idist = {}; % experimental
    for buc = 1:length(rtime)
        rlog2_Idist{1,buc} = log2n_H2BGFP_Int_all{1,buc}-median(log2n_H2BGFP_Int_all{1,1},1);
    end
    
    % Plot final distribution of intensities:
    xrange = [-14:0.25:2];
    for bb = 1:4 %for each time point
        [counts,centers] = hist(clog2_Idist(bb,:),xrange); [datacounts,datacenters] = hist(rlog2_Idist{1,bb},xrange);
        subplot(4,4,bb)
        hold on
        plot(datacenters,datacounts./max(datacounts),'k');
        plot(centers,counts./max(counts),'Color','r','LineStyle',mystyle{sj}); 
        hold off
        xlim([-14 2]); ylim([0 1.1])
        for aja = -1:13; line([-aja -aja],[0 1.1],'Color','k','LineStyle',':'); end
        set(gca,'YTick',[0:0.2:1]); set(gca,'XTick',[-10 -5 0])
        box on
        ylabel('Frequency'); xlabel('log2(H2BGFP int.)')
    end
    
end
