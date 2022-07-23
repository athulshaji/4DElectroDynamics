function MyTestFunction4
clc;clear all;close all;
format long;
Tf=1000
T0=0;
R1=131;% Ohms
C1=0.0107e-9;% Farads
% C1=1e-9;
t=linspace(0,1000,1001);
t=1e-9*t;% converting time into nano-seconds
Pf=ones(1,501)*3.63;% Watts
Pf=[Pf zeros(1,500)];
% Pf=3.63;% Watts
rhs=Pf.*R1.*(1-exp(-t/(R1*C1)))+T0;
figure,plot(t,rhs,'LineWidth',3);grid on;axis tight;

end