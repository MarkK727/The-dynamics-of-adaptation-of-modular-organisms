% Description: This is a comprehensive function 

% Inputs: repeat, regime, FixedParam
% Output: Final

function [simData] = cFGM_simulation( repeat, regime, FixedParam )
%% Set structure "Final" that stores all final meta-information

% Store population mean and variance of (i) phenotypes and (ii) fitness per
% generation for all simulations
simData.DataTable_Averages = cell(1);
simData.DataTable_Variances = cell(1);
simData.Initial_Angle = zeros();
simData.DataTable_Angle = cell(1);

simData.Summary = zeros();

%% Setup key parameters
% N: Population size.
% U: Base rate of mutation per generation.
% R: Rate of Recombination ranging [0, 0.5]

N = regime.N;
U = regime.U;
R = regime.R;

%% << Setup Parameters and variables >>

InitialAngle = FixedParam.InitialAngle;

for i_pos = 1:length(InitialAngle)
% sigW: The standard deviation of the bivariate Gaussian fitness landscape
% We define sigW in terms of multiples of the mutation step size, k. By
% default, we use 10*k to ensure multiple adaptive steps to an optimum.
k = FixedParam.k;
sigW = 10*k;

% Due to unequal selection pressures on modules, the fitness isocline is
% elliptical; here, we find the distance, d, which with an input angle,
% theta, yields the prescribed value of fitness, Default_fitness. By
% default, we start with an intial fitness equivalent to 2 stdv away from
% the optimum.
theta = InitialAngle(i_pos); 
Default = -2*[sigW; sigW];
Default_fitness = exp(-(vecnorm(Default))^2/2/sigW^2);
d = Find_InitialPoint( theta, FixedParam.SelectionBias, Default_fitness, sigW );
d = d(d>0);

WT = d.*[-cos(theta); -sin(theta)];

% [a, b]: Unequal selection pressures on modules
% Complexity: The sub-dimention of each module [n1, n2];
a = FixedParam.SelectionBias(1);
b = FixedParam.SelectionBias(2);
Complexity = FixedParam.Complexity;

%% Repeat 

for i_repeat = 1:repeat
% Population matrix storing meta-information of all individuals.
Popmat = zeros(3,N);

Popmat(1,:) = WT(1); 
Popmat(2,:) = WT(2); 

F = sqrt(a*Popmat(1,:).^2 + b*Popmat(2,:).^2); 
Popmat(3,:) = exp( -F.^2/2/sigW^2 );

%% << Main Loop >>

End = false;

t = 1;
% Ends the simulation at this fitness level.
Final_W  = 0.95;

% Stores average phenotype.
Averages = zeros();
Variances = zeros();

while End == false
    %% Mutation
    
    nMut = poissrnd(N*U); % The number of mutation events
    MutID = datasample(1:N, nMut, 'Replace', false);
    Mutmat = Popmat(:,MutID);

    % Update probability of mutation targeting module X, p
    p = Complexity(1)/sum(Complexity);    
    Target_X = binornd(1, p*ones(1,nMut));
    Target_Y = ~Target_X;
    
    % Calculate the "fitness" of module X, W(d) = exp(-d.^2/2) 
    % given that abs(X) = -log(W(d)).
    wX = exp(-abs(Mutmat(1,:)));
    VAR_X = sqrt( ( 2*(k^2)*(-log(wX))/Complexity(1) ) );
    % N(mu, sig^2) = N(0,1)*sig + mu. DFE.
    Displacement_X = randn(1,nMut).*sqrt(VAR_X) + (-k^2)/2; 
    
    % Calculate the "fitness" of module Y, W(d) = exp(-d.^2/2) 
    % given that abs(Y) = -log(W(d)).
    wY = exp(-abs(Mutmat(2,:)));
    VAR_Y = sqrt( ( 2*(k^2)*(-log(wY))/Complexity(2) ) );
    % N(mu, sig^2) = N(0,1)*sig + mu
    Displacement_Y = randn(1,nMut).*sqrt(VAR_Y) + (-k^2)/2; 
    
    Mutmat(1,:) = log(wX) + Target_X.*Displacement_X;
    Mutmat(2,:) = log(wY) + Target_Y.*Displacement_Y;
    F = sqrt(a*Mutmat(1,:).^2 + b*Mutmat(2,:).^2); 
    Mutmat(3,:) = exp( -F.^2/2/sigW^2 ); % Fitness 
    
    Popmat(:,MutID) = Mutmat;
    
    %% Recombination

    if R ~= 0
        if R == 1  % Full recombination
            nPair = N/2;
        else
            nPair = mybinornd(N,R/2); % The number of recombining pair.
            if 2*nPair > N 
                nPair = N/2;
            end
        end
        
        RecID = randperm(N,2*nPair); % Index of recombining individual.
        Recmat = Popmat(1:2, RecID);
        % Recombine modules. A module's mutation rate does not change. 
        Recmat(2, 1:nPair) = Popmat(2, nPair+1:2*nPair);
        Recmat(2, nPair+1:2*nPair) = Popmat(2, 1:nPair);
        
        % Recalculate fitness
        F = sqrt(a*Recmat(1,:).^2 + b*Recmat(2,:).^2); 
        Recmat(3,:) = exp( -F.^2/2/sigW^2 ); % Fitness 
        
        Popmat(:,RecID) = Recmat;
    end
    
    %% Selection
    
    PROB = Popmat(3,:)/sum(Popmat(3,:));
    
    Smat = zeros(2, N);
    Smat(1,1:N) = 1:N;
    Smat(2,:) = mnrnd(N, PROB);
    
    % Introduce an Indicator variable.
    Indict = 0; 
    for i = 1:max(Smat(2,:))
        Offsprings = Smat(:, Smat(2,:) == i);
    
        Template = Popmat(:, Offsprings(1,:));
        Template = repmat(Template, 1, i);
        
        Popmat(:, Indict+1 : Indict+length(Template(1,:))) = Template;
        
        Indict = Indict+length(Template(1,:));
    end
    
    %% Meta-information
    Averages(1,t) = mean(Popmat(1,:));
    Averages(2,t) = mean(Popmat(2,:));   
    Averages(3,t) = mean(log(Popmat(3,:)));
    
    Variances(1,t) = var(Popmat(1,:));
    Variances(2,t) = var(Popmat(2,:));
    Variances(3,t) = var(log(Popmat(3,:)));
    
    % Stop the loop when mean fitness is Final_W.  
    if  Averages(3,t) >= log(Final_W)
        End = true; 
        
        X_coor = mean(Popmat(1,:));
        Y_coor = mean(Popmat(2,:));  
     
        Angle = atan(abs(Y_coor/X_coor));
        
        Final_t = t;
    end
    
%% Update Time

    t = t + 1;
    
end

Strain.RPV(1, i_repeat) = Angle;
dr = log(Final_W) - log(Default_fitness);
Strain.drdt(1, i_repeat) = dr/Final_t;

simData.DataTable_Averages{i_pos,i_repeat} = Averages(1:3,:);
simData.DataTable_Variances{i_pos,i_repeat} = Variances(1:3,:);

end

Angle_mean = mean(Strain.RPV);
Angle_std = sqrt(var(Strain.RPV));

drdt_mean = mean(Strain.drdt);
drdt_std = sqrt(var(Strain.drdt));

%% Data recording

    simData.DataTable_Angle(1,i_pos) = {Strain.RPV};
    simData.Initial_Angle(1,i_pos) = theta;
    
    simData.Summary(1,i_pos) = Angle_mean;
    simData.Summary(2,i_pos) = Angle_std;
    
    simData.Summary(3,i_pos) = drdt_mean;
    simData.Summary(4,i_pos) = drdt_std;
        
end