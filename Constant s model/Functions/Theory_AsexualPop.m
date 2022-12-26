
function [Theory] = Theory_AsexualPop( regime, FixedParam )
%% Set structure "Final" that stores all final meta-information

Theory.Record = cell(1);

U = regime.U;

% Three initial distance from optimum
ID = [1,2,3];

for i_pos = 1:3
%% << Setup Parameters and variables >>
% Use the same input parameters we used in the 'Simulation' function.
k = FixedParam.k;
sigW = 10*k;

T1 = ID(i_pos);
T2 = flip(ID);
T2 = T2(i_pos);
FixedParam.WT = [-T1*sigW; -T2*sigW];

% Phenotype of the initial population.
WT = FixedParam.WT;
rho = abs(WT(1)) + abs(WT(2));

iter = true;
Record = zeros(2,1);

t = 1;

while iter == true
Record(1:2,t) = WT;
    
X = WT(1)+k;
Y = WT(2)+k;

UX = (abs(X)/rho)*U;
UY = (abs(Y)/rho)*U;

DirectionVec = ([UX;UY]/norm([UX;UY]))*k;

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