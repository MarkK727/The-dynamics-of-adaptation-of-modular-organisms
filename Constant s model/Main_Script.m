% Elapsed time is 1446 seconds. 2874 for 10^5
tic

% Declare rates of recombination
R_List = [1,0]; 
for i_n = 1:length(R_List)

% Number of iteration
repeat = 2;

FixedParam.k = 0.1;
% Set selection coefficients and initial fitness
FixedParam.ds = 0.3; 
FixedParam.Initial_Fitness = 1;

% Parameters defining evolutionary regimes
regime.N = 10^3; % Here, we fix N.
regime.U = 10^-2;
regime.R = R_List(i_n);

fname = sprintf('Data_%d.mat', i_n);
save(fname, 'regime', 'repeat', 'FixedParam');

[simData] = Simulation( repeat, regime, FixedParam );
save(fname, 'simData', '-append');

if regime.R == 1
    [Theory] = Theory_SexualPop( regime, FixedParam );
    save(fname, 'Theory', '-append');
end

if regime.R == 0
    [Theory] = Theory_AsexualPop( regime, FixedParam );
    save(fname, 'Theory', '-append');
end

clearvars -except n i_n fname FixedParam sigW List R_List

end

toc