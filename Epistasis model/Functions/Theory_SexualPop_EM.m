% Description: This is a comprehensive function 

% Inputs:
% Output: 

function [Theory] = Theory_SexualPop_EM( regime, FixedParam )
%% Set structure "Final" that stores all final meta-information

Theory.Record = cell(1);

%% Setup key parameters
% N: Population size.
% U: Base rate of mutation per generation.
% R: Rate of Recombination ranging [0, 0.5]

N = regime.N;
U = regime.U;

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
a = FixedParam.SelectionBias(1);
b = FixedParam.SelectionBias(2);
rho = abs(Default(1)) + abs(Default(2));

iter = true;
Record = zeros(2,1);

t = 1;

while iter == true
Record(1:2,t) = WT;

X = WT(1)+k;
Y = WT(2)+k;

F = sqrt(a*WT(1).^2 + b*WT(2).^2); 
F_X = sqrt(a*(X^2) + b*WT(2).^2); 
F_Y = sqrt(a*WT(1).^2 + b*(Y^2)); 

sX = log(exp( -F_X.^2/2/sigW^2 )/exp( -F.^2/2/sigW^2 )); 
Ux = (abs(X)/rho)*U;
vX = (sX^2)*( (2*log(N*sX) - log(sX/Ux))/((log(sX/Ux))^2) );
dXdt = (vX/sX)*k;

sY = log(exp( -F_Y.^2/2/sigW^2 )/exp( -F.^2/2/sigW^2 ));    
Uy = (abs(Y)/rho)*U;
vY = (sY^2)*( (2*log(N*sY) - log(sY/Uy))/((log(sY/Uy))^2) );
dYdt = (vY/sY)*k;

DirectionVec = [dXdt;dYdt];

WT = WT + DirectionVec;

if exp( -F.^2/2/sigW^2 ) >= 0.95
   iter = false;
end
t = t+1;
end

%% Data recording

Theory.Record{i_pos, 1} = Record;

end
        
end