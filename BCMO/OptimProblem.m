%% Generalized Rastrigin Function
function f = OptimProblem(x)
%f=(sum(x.^2-10*cos(x.*2*pi)+10)); % Limited Bounds: [-5.12 5.12]; Min = 0
f=(sum(-x.*sin(sqrt(abs(x))))); % Limited Bounds: [-500 500]; Min = ?418.9829+e5
end   