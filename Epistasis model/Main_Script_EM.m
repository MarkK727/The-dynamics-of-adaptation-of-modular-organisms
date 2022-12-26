% This is a comprehensive script which generates all simulation results for
% epistasis model

tic
% Declare rates of recombination
R_List = [1,10^-2,0]; 

for i_n = 1:3

% Number of iteration
repeat = 250;

FixedParam.k = 0.1;
FixedParam.SelectionBias = [1, 2]; 

% Set the initial phenotypic positions described in terms of angles
% pi/24 ~ 11pi/24 rad
% Remove 0 and 90 degrees (on the axes)
InitialAngle = linspace(0,pi/2,13); 
InitialAngle = InitialAngle(2:length(InitialAngle)-1); 

FixedParam.InitialAngle = InitialAngle;

% Parameters defining evolutionary regimes
regime.N = 10^4; % Here, we fix N.
regime.U = 10^-3;
regime.R = R_List(i_n);

fname = sprintf('Data_%d.mat', i_n);
save(fname, 'regime', 'repeat', 'FixedParam');

[simData] = EM_simulation( repeat, regime, FixedParam );
save(fname, 'simData', '-append');

if regime.R == 0
    [Theory] = Theory_AsexualPop_EM( regime, FixedParam );
    save(fname, 'Theory', '-append');
end

if regime.R == 1
    [Theory] = Theory_SexualPop_EM( regime, FixedParam );
    save(fname, 'Theory', '-append');
end

clearvars -except n i_n fname FixedParam sigW List R_List

end

toc