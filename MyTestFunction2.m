
function MyTestFunction2
clc;close all;clear variables;
format long
global Vd Is Rs eta Vth
Vd=0.6;Is=2.5e-9;Rs=0.6;eta=1.7;Vth=0.026;
x0=0;
options = optimoptions('fsolve',...
    'Display','iter-detailed',...
    'Algorithm','levenberg-marquardt',...
    'FunctionTolerance',1e-7,...
    'MaxIterations',1000,...
    'MaxFunctionEvaluations',100000);
[Id,fval,exitflag,output]=fsolve(@fun,x0,options);
[Id (exitflag)]
end

function error=fun(Id)
global Is Vd Rs eta Vth
error=Is*(exp((Vd-Rs*Id)/(eta*Vth))-1)-Id;
end

