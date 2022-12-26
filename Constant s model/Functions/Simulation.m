function [simData] = Simulation( repeat, regime, FixedParam )

% Store population mean and variance of (i) phenotypes and (ii) fitness per
% generation for all simulations
simData.DataTable_Averages = cell(1);

%% Setup key parameters
% N: Population size.
% U: Base rate of mutation per generation.
% R: The fraction of pairs in the populution going through recombination.
% Numerically, it is the cannonical rate of recombination [0, 0.5]
% mutliplied by 2.

N = regime.N;
U = regime.U;
R = regime.R;

% Three initial distance from optimum
ID = [1,2,3];

%%
for i_pos = 1:3
%% << Setup Parameters and variables >>
% sigW: The standard deviation of the bivariate Gaussian fitness landscape
% We define sigW in terms of multiples of the mutation step size, k. By
% default, we use 10*k to ensure multiple adaptive steps to an optimum.
k = FixedParam.k;
sigW = 10*k;
    
% Set strength of selection
ds = FixedParam.ds;
Initial_Fitness = FixedParam.Initial_Fitness;

% Initial Population Phenotype (Trait)
T1 = ID(i_pos);
T2 = flip(ID);
T2 = T2(i_pos);
FixedParam.WT = [-T1*sigW; -T2*sigW];
    
%% << Create transition matrix >>
% We simulate mutation using transition probability matrix containing
% beneficial mutation rates (UMat) at each phenotypic coordinate. The final
% vectors, p1vec and p2vec, is transition probability matrix for module x
% and module y. 

% Phenotype of the initial founder.
WT = FixedParam.WT;

% Assign finite mutation supply
Xn = abs(WT(1))/k; % Total mutation supply of trait X
Yn = abs(WT(2))/k; % Total mutation supply of trait Y

MP_X = zeros(Xn+1,2); % Mutation Profile of Trait X
MP_X(:,1) = 0:Xn; % Fisrt row represents the number of mutation on the trait X
MP_X(1,2) = N; % The wild type has no mutation. 

MP_Y = zeros(Yn+1,2); % Mutation Profile of Trait Y
MP_Y(:,1) = 0:Yn; % Fisrt row represents the number of mutation on the trait Y
MP_Y(1,2) = N; % The wild type has no mutation. 

dX = MP_X(:,1).*k; % Phenotypic displacement of X
dY = MP_Y(:,1).*k; % Phenotypic displacement of Y

FGM_X = WT(1) + dX';
FGM_Y = WT(2) + dY;

Popmat = zeros(Yn+1, Xn+1); % Number of individuals in strain (i,j)
[row,col] = size(Popmat);

X = repmat((MP_X(:,1).*ds)', row, 1);
Y = repmat((MP_Y(:,1).*ds), 1, col);
FitnessMat = X + Y;

FitnessMat = exp(FitnessMat);
Fitnessvec = FitnessMat(:);

UX = repmat(Xn-(MP_X(:,1))', row, 1);
UY = repmat(Yn-(MP_Y(:,1)), 1, col);
UMat = ((UX)./(Xn + Yn))*U + ((UY)./(Xn + Yn))*U;
Uvec = UMat(:);

p1 = ( (UX./(Xn + Yn))*U )./UMat;
p2 = ( (UY./(Xn + Yn))*U )./UMat;
p1vec = p1(:);
p2vec = p2(:);

%% Main Loop

for i_repeat = 1:repeat
clear Popmat

Popmat = zeros(Yn+1, Xn+1); % Number of individuals in strain (i,j)
Popmat(1,1) = N;
Popvec = Popmat(:); % Change to a single vector for computational purpose

%% << Loop Variable >>

End = false;

End_F = log(Initial_Fitness) + (Xn-2)*ds + (Yn-2)*ds;
% Ends the simulation at this fitness level.
Final_W  = exp(End_F); 

% Reset time
t = 1;

% Stores average phenotype.
Averages = zeros();

while End == false    
    %% Mutation
    
    Mutvec = mybinornd(Popvec, Uvec);
    dX_Pop = mybinornd(Mutvec, p1vec);
    dY_Pop = mybinornd(Mutvec, p2vec);
    
    Popvec = Popvec - Mutvec;
    Popmat = reshape(Popvec, [row,col]);
    
    dX_Pop = reshape(dX_Pop, [row,col]);
    dX_Pop(:,2:col) = dX_Pop(:,1:col-1);
    dX_Pop(:,1) = 0;
    
    dY_Pop = reshape(dY_Pop, [row,col]);
    dY_Pop(2:row,:) = dY_Pop(1:row-1,:);
    dY_Pop(1,:) = 0;
    
    Popvec = Popvec + dX_Pop(:) + dY_Pop(:);
    
    %% Recombination
    % The effect of free recombination between unlinked Xloci and Yloci are
    % approximated implicitly, by re-drawing frequencies of phenotypes
    % based on the current frequencies of (X, Y).
   
    if R ~= 0
        Popmat = reshape(Popvec, [row,col]);
        f_X = (sum(Popmat)/sum(Popmat(:))); % Frequencies of X phenotypes
        f_Y = (sum(Popmat,2)/sum(Popmat(:))); % Frequencies of Y phenotypes
        f_XY = f_X.*f_Y;

        RecombPop = poissrnd(N*f_XY(:));
        RecombPop = reshape(RecombPop, [row,col]);
        Popvec = RecombPop(:);
    end
    
    %% Selection
    
    PROB = (Popvec.*Fitnessvec)/sum(Popvec.*Fitnessvec);
    meanfitness = sum(Popvec.*Fitnessvec)/sum(Popvec);
    
    Popvec = poissrnd(PROB*N); % Poission estimation of multinomial draws.
    
    %% Meta-information
    
    X_coor = sum(sum(Popmat).*FGM_X)/sum(Popmat(:));
    Y_coor = sum(sum(Popmat,2).*FGM_Y)/sum(Popmat(:));
    
    Averages(1,t) = X_coor;
    Averages(2,t) = Y_coor;
    Averages(3,t) = log(meanfitness);
    
    if  Averages(3,t) >= log(Final_W)
        End = true;
    end
    
%% Update Time

    t = t + 1;
    
end

simData.DataTable_Averages{i_pos,i_repeat} = Averages(1:3,:);

end
    
end