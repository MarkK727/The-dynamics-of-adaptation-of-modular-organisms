function d = Find_InitialPoint( theta, ShapeParameter, Default_fitness, sigW )

if nargin < 3
    Default_fitness = 0.000123409804086679;
    sigW = 1;
end

a1 = ShapeParameter(1);
a2 = ShapeParameter(2); 

syms r
eqn = Default_fitness == exp(-sqrt(a1*(r*cos(theta))^2 + a2*(r*sin(theta))^2)^2/2/sigW^2); 
sol = solve(eqn, r);
d = double(sol);

