% Description: This is a comprehensive function 

% Inputs:
% Output: 

function [Opt] = OptimalPath( FixedParam )
%% Set structure "Final" that stores all final meta-information

Opt.Path = cell(1);

%% << Setup Parameters and variables >>

InitialAngle = FixedParam.InitialAngle;

% [a, b]: Unequal selection pressures on modules
% rho: Scaling factor for the basal mutation rate
a = FixedParam.SelectionBias(1);
b = FixedParam.SelectionBias(2);

% Take gradient function
syms x1 x2
f(x1,x2) = exp(-sqrt(FixedParam.SelectionBias(1)*x1^2 + FixedParam.SelectionBias(2)*x2^2)^2/2/1^2);
grad = gradient(f,[x1,x2]);

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

iter = true;
Path = zeros(2,1);

t = 1;

while iter == true
Path(1:2,t) = WT;

gVec = grad(WT(1),WT(2));
gVec = double(gVec);
gVec = (gVec/vecnorm(gVec))*k;

F = sqrt(a*WT(1).^2 + b*WT(2).^2); 

WT = WT + gVec;

if exp( -F.^2/2/sigW^2 ) >= 0.95
   iter = false;
end
t = t+1;
end

%% Data recording

Opt.Path{i_pos, 1} = Path;

end
        
end