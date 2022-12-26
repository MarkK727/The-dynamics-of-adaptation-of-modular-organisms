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
FixedParam.Complexity = [20;20];

% Set the initial phenotypic positions described in terms of angles
% pi/24 ~ 11pi/24 rad
% Remove 0 and 90 degrees
InitialAngle = linspace(0,pi/2,13); 
InitialAngle = InitialAngle(2:length(InitialAngle)-1); 

FixedParam.InitialAngle = InitialAngle;

% Parameters defining evolutionary regimes
regime.N = 10^4; % Here, we fix N.
regime.U = 10^-3;
regime.R = R_List(i_n);

fname = sprintf('Data_%d.mat', i_n);
save(fname, 'regime', 'repeat', 'FixedParam');

[simData] = cFGM_simulation( repeat, regime, FixedParam );
save(fname, 'simData', '-append');

clearvars -except n i_n fname FixedParam sigW List R_List

end

%% Aesthetics

figure(1)

xticks(0*pi:(1/12)*pi:(1/2)*pi)
xticklabels({'0\pi','\pi_{/12}','\pi_{/6}','\pi_{/4}','\pi_{/3}', '5\pi_{/12}', '\pi_{/2}'});

yticks(0*pi:(1/12)*pi:(1/2)*pi)
yticklabels({'','\pi_{/12}','\pi_{/6}','\pi_{/4}','\pi_{/3}', '5\pi_{/12}', '\pi_{/2}'});

H=gca;
H.LineWidth=0.8; %change to the desired value 

xlabel('Initial Angles')
ylabel('Final Angles')

set(gca, 'FontName', 'Times New Roman', 'FontSize',15)

figure(2)

xticks(0*pi:(1/12)*pi:(1/2)*pi)
xticklabels({'0\pi','\pi_{/12}','\pi_{/6}','\pi_{/4}','\pi_{/3}', '5\pi_{/12}', '\pi_{/2}'});

H=gca;
H.LineWidth=0.8; %change to the desired value 

xlabel('Initial Angles')
ylabel('Average rate of evolution (dr/dt)')

set(gca, 'FontName', 'Times New Roman', 'FontSize',15)

toc