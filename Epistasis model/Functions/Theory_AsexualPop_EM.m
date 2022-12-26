% Description: This is a comprehensive function 

% Inputs:
% Output: 

function [Theory] = Theory_AsexualPop_EM( regime, FixedParam )
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
a = FixedParam.SelectionBias(1);
b = FixedParam.SelectionBias(2);
rho = abs(Default(1)) + abs(Default(2));

iter = true;
Record = zeros(2,1);

Root1 = zeros(1,1);
Root2 = zeros(1,1);

t = 1;


while iter == true
Record(1:2,t) = WT;
    
X = WT(1)+k;
Y = WT(2)+k;

d0 = sqrt(a*WT(1).^2 + b*WT(2).^2);
dX = sqrt(a*X.^2 + b*WT(2).^2);
dY = sqrt(a*WT(1).^2 + b*Y.^2);

W0 = exp( -d0.^2/2/sigW^2 );
WX = exp( -dX.^2/2/sigW^2 );
WY = exp( -dY.^2/2/sigW^2 );

sX = log(WX/W0);
sY = log(WY/W0);

UX = (abs(X)/rho)*U;
UY = (abs(Y)/rho)*U;

syms s0 U0 N0
AdaptRate = ( (s0^2)*( (2*log(N0*s0))-log(s0/U0) ) )/( log(s0/U0) )^2;

%vX = double(subs(AdaptRate, {s0, U0, N0}, {sX, UX, N}));
vY = double(subs(AdaptRate, {s0, U0, N0}, {sY, UY, N}));

eqn = subs(AdaptRate, {s0, N0}, {sX, N}) == vY;
UY_tilde = double(solve(eqn,U0));

Root1(t) = UY_tilde(1);
Root2(t) = UY_tilde(2);

UY_tilde = UY_tilde(sX./UY_tilde > 1); % s/Ub >> 1

DirectionVec = ([1;UY_tilde/UX]/norm([1;UY_tilde/UX]))*k;

WT = WT + DirectionVec;

F = sqrt(a*WT(1).^2 + b*WT(2).^2); 

if exp( -F.^2/2/sigW^2 ) >= 0.95
   iter = false;
end

t = t+1;

end

%% Data recording

Theory.Record{i_pos, 1} = Record;
Theory.Record{i_pos, 2} = Root1;
Theory.Record{i_pos, 3} = Root2;

end
        
end