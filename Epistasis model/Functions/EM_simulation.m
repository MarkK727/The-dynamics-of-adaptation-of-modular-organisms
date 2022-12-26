% Description: This is a comprehensive function 

% Inputs: repeat, regime, FixedParam
% Output: Final

function [simData] = EM_simulation( repeat, regime, FixedParam )
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
% theta, yields the prescribed value of fitness, Default_fitness.
theta = InitialAngle(i_pos); 
Default = -2*[sigW; sigW];
Default_fitness = exp(-(vecnorm(Default))^2/2/sigW^2);
d = Find_InitialPoint( theta, FixedParam.SelectionBias, Default_fitness, sigW );
d = d(d>0);

WT = d.*[-cos(theta); -sin(theta)];

% [a, b]: Unequal selection pressures on modules
% rho: Scaling factor for the basal mutation rate
rho = abs(Default(1)) + abs(Default(2));
a = FixedParam.SelectionBias(1);
b = FixedParam.SelectionBias(2);

%% Repeat 

for i_repeat = 1:repeat
% Population matrix storing meta-information of all individuals.
Popmat = zeros(5,N);

Popmat(1,:) = WT(1); % Module X
Popmat(2,:) = WT(2); % Module Y

F = sqrt(a*Popmat(1,:).^2 + b*Popmat(2,:).^2); 
Popmat(3,:) = exp( -F.^2/2/sigW^2 ); % Fitness 

rho_X = abs(Popmat(1,:))/rho;
rho_Y = abs(Popmat(2,:))/rho;
Popmat(4,:) = rho_X.*U; % Mutation rate
Popmat(5,:) = rho_Y.*U; % Mutation rate

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
    
    nMut_X = poissrnd(sum(Popmat(4,:))); % The number of mutation events
    MutID_X = datasample(1:N, nMut_X, 'Weights', Popmat(4,:),'Replace',false);
    Mutmat_X = Popmat(:,MutID_X);
    Mutmat_X(1,:) = Mutmat_X(1,:) + k;
    
    nMut_Y = poissrnd(sum(Popmat(5,:))); % The number of mutation events
    MutID_Y = datasample(1:N, nMut_Y, 'Weights', Popmat(5,:),'Replace',false);
    Mutmat_Y = Popmat(:,MutID_Y);
    Mutmat_Y(2,:) = Mutmat_Y(2,:) + k;
    
    F = sqrt(a*[Mutmat_X(1,:),Mutmat_Y(1,:)].^2 + b*[Mutmat_X(2,:),Mutmat_Y(2,:)].^2); 
    
    Mutmat_X(3,:) = exp( -F(1:length(Mutmat_X(3,:))).^2/2/sigW^2 ); 
    Mutmat_Y(3,:) = exp( -F(length(Mutmat_X(3,:))+1:end).^2/2/sigW^2 ); 

    Mutmat_X(4,:) = (abs(Mutmat_X(1,:))/rho).*U;
    Mutmat_Y(5,:) = (abs(Mutmat_Y(2,:))/rho).*U;
    
    Popmat(:,MutID_X) = Mutmat_X;
    Popmat(:,MutID_Y) = Mutmat_Y;
    
    %% Recombination

    if R ~= 0
        if R == 1  % Full recombination
            nPair = N/2;
        else
            nPair = mybinornd(N, R/2); % The number of recombining pair.
            if 2*nPair > N 
                nPair = N/2;
            end
        end
        
        RecID = randperm(N,2*nPair); % Index of recombining individual.
        Recmat = Popmat(1:2, RecID);
        Recmat(2, 1:nPair) = Popmat(2, nPair+1:2*nPair);
        Recmat(2, nPair+1:2*nPair) = Popmat(2, 1:nPair);
        
        F = sqrt(a*Recmat(1,:).^2 + b*Recmat(2,:).^2); 
        Recmat(3,:) = exp( -F.^2/2/sigW^2 ); % Fitness 
        Recmat(4,:) = (abs(Recmat(1,:))/rho).*U;
        Recmat(5,:) = (abs(Recmat(2,:))/rho).*U;
        
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