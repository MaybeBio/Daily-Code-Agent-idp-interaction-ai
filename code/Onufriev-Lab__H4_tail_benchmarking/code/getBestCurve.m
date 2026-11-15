
function best_curve = getBestCurve(dir,nRuns,curvesPerEnsemble)
%% Take output from GAJOE v2-3 and extract PDB numbers in all ensembles
%   Essentially recreates the best_curve.txt file generated from gajoe 1.3
%   that was removed from gajoe 2 for SOME reason.
%
%   First run GAJOE v2-3. See https://www.embl-hamburg.de/biosaxs/manuals/eom.html   
%   Note: you need to select yes to the option:
%   > Save workfiles with the entire GA runs? (default: no) ................ : y
%
%       GW - May 2024 
%


%% 
GAresult = {};
selectedStructures = [];

for i = 1:nRuns
    raw = readmatrix([dir,'/GA_CYCLE_',num2str(i),'.txt'],"Delimiter",'-');
    rawActual = raw(end-(curvesPerEnsemble-1):end,:);
    rawStructures = reshape(rawActual,[],1);
    rawStructures = rawStructures(~isnan(rawStructures));

    selectedStructures = [selectedStructures; rawStructures];
    
end

selectedStructures = categorical(selectedStructures);
table = tabulate(selectedStructures);
for j = 1:numel(table(:,1))
    PDBnums(j) = str2double(table(j,1));
    justGetTheFuckingValue = table(j,2);
    freqnums(j)= justGetTheFuckingValue{1};
end
PDBnums = PDBnums';
freqnums = freqnums';

best_curve = [PDBnums , freqnums];


end
