function MyFunctionSPICE
% written on 11-Apr-2022 to verify the SPICE code
clc;close all;clear;format shortEng;

Cap=1e-9;
Res=500;
freq=5.2e9;
start_time=0; end_time=5e-9;
delt=(1/freq)/100;
Vpre=0;
index=0;
MNA=[1/Res -1/Res 1;...
    -1/Res 1/Res 0;...
     0 Cap/delt -1];
tic;[L,U,P]=lu(MNA);toc;
tic;
for time=start_time:delt:end_time
    RHS=[0;sin(2*pi*freq*time);(Cap/delt)*Vpre];
    V=U\(L\(P*RHS));
    Vpre=V(2);
    index=index+1;
    V1(index)=V(1);V2(index)=V(2);
end
toc;
figure(1),plot(0:delt:end_time,V1);grid on;figure(2),plot(0:delt:end_time,V2);grid on;
end