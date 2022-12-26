
function [Theory] = Theory_SexualPop( regime, FixedParam )
%% Set structure "Final" that stores all final meta-information

Theory.Record = cell(1);

%% << Setup Parameters and variables >>
% Use the same input parameters we used in the 'Simulation' function.

N = regime.N;
U = regime.U;

% Three initial distance from optimum
ID = [1,2,3];

for i_pos = 1:3

k = FixedParam.k;
sigW = 10*k;

T1 = ID(i_pos);
T2 = flip(ID);
T2 = T2(i_pos);
FixedParam.WT = [-T1*sigW; -T2*sigW];

% Phenotype of the initial founder.
WT = FixedParam.WT;
% rho: Scaling factor for the basal mutation rate
rho = abs(WT(1)) + abs(WT(2));

iter = true;
Record = zeros(2,1);

t = 1;

while iter == true

Record(1:2,t) = WT;

X = WT(1)+k;
Y = WT(2)+k;

sX = FixedParam.ds;

Ux = (abs(X)/rho)*U;
vX = (sX^2)*( (2*log(N*sX) - log(sX/Ux))/((log(sX/Ux))^2) );
dXdt = (vX/sX)*k;

sY = FixedParam.ds;

Uy = (abs(Y)/rho)*U;
vY = (sY^2)*( (2*log(N*sY) - log(sY/Uy))/((log(sY/Uy))^2) );
dYdt = (vY/sY)*k;

DirectionVec = [dXdt;dYdt];

WT = WT + DirectionVec;

if WT(1) > -1.1*k || WT(2) > -1.1*k
   iter = false;
end
t = t+1;
end

%% Data recording

Theory.Record{i_pos, 1} = Record;

end
        
end