
%% putting together GAJOE results
%first, get end-to-end distance R by running pdb_Rmeas.m in directory with all the PDBs of full pool, and 
% load in the Ree values into the workspace
%
% GW Feb 2023 - removed EFRET plotting and made plots publication quality 
% GW May 2024 - adapted for use with GAJOE v2.1, which does not output a
% best_curve txt file anymore. 
% 
%
clear; close all

%% User-defined parameters
dir = 'ExampleData/H4_tail_WT_ff19SB_OPC/GA001/curve_1';
nRuns = 100; % What you entered for: Number of times you want the genetic algorithm repeated (min. 1, default: 100)
curvesPerEnsemble = 20; % What you entered for: Maximum number of curves per ensemble (min. 1, max. 50) 

xlim_Dmax = [20 80]; ylim_Rg = [5 30]; xlim_R = [0 80]; % Plot axes limits 

% Plot histogram bins and file sizes 
gajoeHistBins = 15;
fontSize = 33;
axisFontSize = 22;
binsRg = 22; 
binsDmax = 22; 
binsR = 22;
% -- -- -- -- -- -- -- -- -- -- -- --


% Load best_curve 
bestcurve = getBestCurve(dir,nRuns,curvesPerEnsemble);

%% Import variables 
% R (Ree), Rg, Dmax, compute E (fret eff) from R
load([dir,'/Ree_values.mat']) % Loads variable 'R'
Sizelist = importdata([dir,'/Size_list.txt']);
pool.dmax = Sizelist(:,3);
pool.rg = Sizelist(:,2);
pool.r = R(1:end);

% GAJOE results
gaj = length(bestcurve(:,1));


%% Delineate GAJOE structures
Dmax_gaj = []; Rg_gaj = []; R_gaj = []; E_gaj = [];
for n = 1:gaj
    D_Rep = repmat(Sizelist(bestcurve(n,1),3),bestcurve(n,2),1); 
    Rg_Rep = repmat(Sizelist(bestcurve(n,1),2),bestcurve(n,2),1);
    R_Rep = repmat(R(bestcurve(n,1)),bestcurve(n,2),1);
    Dmax_gaj = [Dmax_gaj; D_Rep]; Rg_gaj = [Rg_gaj; Rg_Rep];
    R_gaj = [R_gaj; R_Rep];
    clear D_Rep Rg_Rep R_rep 
end
gajoe.dmax = Dmax_gaj; gajoe.rg = Rg_gaj; gajoe.r = R_gaj;


%% Plot 2D-parameterized (Rg, R) structural pool 

FigHandle = figure('Position', [350, 250, 550, 550], 'Name', '2D data');
set(gcf, 'color', 'w')
box on
grid on
set(gcf,'color','w')


subplot(3,3,[4,5,7,8])
hold all
set(gca,'LineWidth',2)
set(gca,'FontSize',20)
set(gca,'FontSize',axisFontSize)
grid on; box on
xlabel('$R_{ee}$ ($\AA$)','FontSize',fontSize,'Interpreter','latex') 
ylabel('$R_{g}$ ($\AA$)','FontSize',fontSize,'Interpreter','latex')
xlim(xlim_R)
ylim(ylim_Rg)

plot(pool.r,pool.rg,'.','Color',[0.8500 0.3250 0.0980]); % full MD pool
scatter(gajoe.r,gajoe.rg,'o','MarkerEdgeColor',[0.3010 0.7450 0.9330],'LineWidth',2); % EOM pool

legend('full pool','selected by EOM','Location','NorthWest')


% Plot R subhistograms of full pool and EOM pool
subplot(3,3,[1,2])
xlim(xlim_R)
grid on; box on
set(gca,'LineWidth',1.5)
set(gca,'FontSize',20)
hold all
hpoolr = histogram(pool.r, binsR,'FaceColor',[0.8500 0.3250 0.0980]);
hpoolr.BinCounts = hpoolr.BinCounts/max(hpoolr.BinCounts); % GW added Sept 2025 - normalize histograms

hgajr = histogram(gajoe.r, gajoeHistBins, 'FaceColor',[0.3010 0.7450 0.9330]);
hgajr.BinCounts = hgajr.BinCounts/max(hgajr.BinCounts); % GW added Sept 2025 - normalize histograms


% Plot Rg subhistograms of full pool and EOM pool
subplot(3,3,[6,9])
set(gca,'xdir','reverse')
xlim(ylim_Rg)
grid on; box on
set(gca,'LineWidth',1.5)
set(gca,'FontSize',20)
hold all
hpool = histogram(pool.rg, binsRg,'FaceColor',[0.8500 0.3250 0.0980]);
hpool.BinCounts = hpool.BinCounts/max(hpool.BinCounts); % GW added Sept 2025 - normalize histograms

hgaj = histogram(gajoe.rg, gajoeHistBins, 'FaceColor',[0.3010 0.7450 0.9330]);
hgaj.BinCounts = hgaj.BinCounts/max(hgaj.BinCounts); % GW added Sept 2025 - normalize histograms
camroll(-90)




%% Print mean values
['Mean Rg of full pool: ',num2str(mean(pool.rg))]
['Standard error in Rg of full pool: ',num2str(std(pool.rg)./sqrt(numel(pool.rg)))]
['Mean Rg of EOM pool: ',num2str(mean(gajoe.rg))]
['Standard error in Rg of EOM pool: ',num2str(std(gajoe.rg)./sqrt(numel(gajoe.rg)))]

% ['Mean Dmax of full pool: ',num2str(mean(pool.dmax))]
% ['Standard error in Dmax of full pool: ',num2str(std(pool.dmax)./sqrt(numel(pool.dmax)))]
% ['Mean Dmax of EOM pool: ',num2str(mean(gajoe.dmax))]
% ['Standard error in Dmax of EOM pool: ',num2str(std(gajoe.dmax)./sqrt(numel(gajoe.dmax)))]

['Mean Ree of full pool: ',num2str(mean(pool.r))]
['Standard error in Ree of full pool: ',num2str(std(pool.r)./sqrt(numel(pool.r)))]
['Mean Ree of EOM pool: ',num2str(mean(gajoe.r))]
['Standard error in Ree of EOM pool: ',num2str(std(gajoe.r)./sqrt(numel(gajoe.r)))]



%% Plot 2D-parameterized (Rg, Dmax) structural pool 

FigHandle = figure('Position', [350, 250, 550, 550], 'Name', '2D data');
set(gcf, 'color', 'w')
box on
grid on
set(gcf,'color','w')


subplot(3,3,[4,5,7,8])
hold all
set(gca,'LineWidth',2)
set(gca,'FontSize',20)
set(gca,'FontSize',axisFontSize)
grid on; box on
xlabel('$D_{max}$ ($\AA$)','FontSize',fontSize,'Interpreter','latex') 
ylabel('$R_{g}$ ($\AA$)','FontSize',fontSize,'Interpreter','latex')
xlim(xlim_Dmax) 
ylim(ylim_Rg)

plot(pool.dmax,pool.rg,'.','Color',[0.8500 0.3250 0.0980]); % full MD pool
scatter(gajoe.dmax,gajoe.rg,'o','MarkerEdgeColor',[0.3010 0.7450 0.9330],'LineWidth',2); % EOM pool

legend('full pool','selected by EOM','Location','NorthWest')


% Plot Dmax subhistograms of full pool and EOM pool
subplot(3,3,[1,2])
xlim(xlim_Dmax)
grid on; box on
set(gca,'LineWidth',1.5)
set(gca,'FontSize',20)
hold all
histogram(pool.dmax, binsDmax,'FaceColor',[0.8500 0.3250 0.0980])
histogram(gajoe.dmax, 'FaceColor',[0.3010 0.7450 0.9330])

% Plot Rg subhistograms of full pool and EOM pool
subplot(3,3,[6,9])
set(gca,'xdir','reverse')
xlim(ylim_Rg)
grid on; box on
set(gca,'LineWidth',1.5)
set(gca,'FontSize',20)
hold all
histogram(pool.rg, binsRg,'FaceColor',[0.8500 0.3250 0.0980])
histogram(gajoe.rg, gajoeHistBins, 'FaceColor',[0.3010 0.7450 0.9330])
camroll(-90)




%% Calculate JS divergence between MD and EOM ensemble distributions 
% GW added Sept 2025
% Method from plotIterativeAnalysis.m (RNA condensates software) 

figure('Visible','off') 
% For Rg 
h_mdpool = histogram(pool.rg, binsRg);
mdpool = histcounts(pool.rg, binsRg, 'Normalization','pdf');
%mdpool = mdpool ./ max(mdpool); % normalization 

h_eompool = histogram(gajoe.rg, binsRg);
eompool = histcounts(gajoe.rg, binsRg, 'Normalization','pdf');
%eompool = eompool ./ max(eompool); % normalization 

p1 = mdpool;
p2 = eompool;

JSD1 = H(0.5.*p1 + 0.5.*p2);
JSD2 = 0.5.*H(p1);
JSD3 = 0.5.*H(p2);

JSD_rg = JSD1 - JSD2 - JSD3


% For R (R_EE)
h_mdpool = histogram(pool.r, binsR);
mdpool = histcounts(pool.r, binsR, 'Normalization','pdf');
%mdpool = mdpool ./ max(mdpool); % normalization 

h_eompool = histogram(gajoe.r, binsR);
eompool = histcounts(gajoe.r, binsR, 'Normalization','pdf');
%eompool = eompool ./ max(eompool); % normalization 

p1 = mdpool;
p2 = eompool;

JSD1 = H(0.5.*p1 + 0.5.*p2);
JSD2 = 0.5.*H(p1);
JSD3 = 0.5.*H(p2);

JSD_r = JSD1 - JSD2 - JSD3



% Subfunction to calculate Shannon entropy 
function H = H(p)
    body = p.*log(p); 
    body(isnan(body)) = 0;
    H = -sum(body);
end

