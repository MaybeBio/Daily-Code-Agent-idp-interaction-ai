
clear; close all

%% Plot EOM fit to SAXS data 
%   GW - February 2023
% Scale factors are simply to visually scale the experimental and fitted curves; they are mapped out in the EOMfitScaling RAW workspace

%% User-defined parameters 
folder = ['ExampleData/H4_tail_WT_ff19SB_OPC/GA001/curve_1'];

% Define plot colour 
theColour = [238 68 0]./255; %(used for H4WT)
%theColour = [46 80 122]./255; %(used for H4KQ)

expData = readmatrix(['SAXSData/','H4WT_150Na.dat'],'NumHeaderLines',4);
fitFile = readmatrix([folder,'/','profiles_001_1.fit'],'FileType','text','NumHeaderLines',1);

% Axis limits for plots 
xlim_SAXS = [0 0.3];
ylim_residuals = [-1E-2 1E-2];

% -- -- -- -- -- -- -- -- -- -- -- --


%% SAXS vs EOM 

%expData = [fitFile(:,1), fitFile(:,2), fitFile(:,3)];
crysolEOM = [fitFile(:,1), fitFile(:,4)];

figure
box on
grid on
set(gcf,'color','w')
set(gca,'LineWidth',2)
set(gca,'FontSize',20)
xlabel('q ($\AA^{-1}$)','Interpreter','latex','FontSize',20)
ylabel('I(q)','Interpreter','latex','FontSize',20)

hold on
errorbar(expData(:,1),expData(:,2),expData(:,3),'Color',[0.7 0.7 0.7])
plot(crysolEOM(:,1),crysolEOM(:,2),'Color',theColour,'LineWidth',3)

xlim(xlim_SAXS)

legend('SAXS experiment','EOM ensemble')

%% SAXS vs EOM , Kratky form 
figure('Units', 'Normalized','Position',[0.3182    0.0852    0.3651    0.8120])
subplot(2,1,1)
box on
grid on
set(gcf,'color','w')
set(gca,'LineWidth',2)
set(gca,'FontSize',20)
xlabel('q ($\AA^{-1}$)','Interpreter','latex','FontSize',20)
ylabel('$q^{2}$ I(q)','Interpreter','latex','FontSize',20)

hold on
errorbar(expData(:,1),expData(:,1).^2.*expData(:,2),expData(:,1).^2.*expData(:,3),'Color',[0.7 0.7 0.7],'Marker','.','LineWidth',1.5)
plot(crysolEOM(:,1),crysolEOM(:,1).^2.*crysolEOM(:,2),'Color',theColour,'LineWidth',3)

xlim(xlim_SAXS)

legend('SAXS experiment','EOM ensemble','Location','SouthEast')

%% Residuals  
subplot(2,1,2)
box on
grid on
set(gcf,'color','w')
set(gca,'LineWidth',2)
set(gca,'FontSize',20)
xlabel('q ($\AA^{-1}$)','Interpreter','latex','FontSize',20)
ylabel('$I_{exp}$ - $I_{theor}$','Interpreter','latex','FontSize',20)

xPoints_crysol = crysolEOM(:,1);
expData_interp_crysol = interp1(expData(:,1),expData(:,2), xPoints_crysol); 
expErrs_interp_crysol = interp1(expData(:,1),expData(:,3), xPoints_crysol);

residualsEOM = expData_interp_crysol - crysolEOM(:,2);

hold all
plot(xPoints_crysol,residualsEOM,'o-','Color',theColour,'LineWidth',2)


expErrs_interp_crysol(isnan(expErrs_interp_crysol))=0;
patch([xPoints_crysol; flipud(xPoints_crysol)], [0+expErrs_interp_crysol; flipud(0-expErrs_interp_crysol)], [0.7 0.7 0.7], 'FaceAlpha',0.5, 'EdgeColor','none')
% crysolErrs = crysolEOM(:,3);
% patch([xPoints_crysol; flipud(xPoints_crysol)], [0+crysolErrs; flipud(0-crysolErrs)], [0.8 0.8 0.8], 'FaceAlpha',0.5, 'EdgeColor','none')


plot(xPoints_crysol,0.*xPoints_crysol,'k--','LineWidth',3)

xlim(xlim_SAXS)
ylim(ylim_residuals)

legend('EOM ensemble','Error margin')





