function MyFunctionSPICE_diode_series_RC_nw
% written on 11-Apr-2022 to verify the SPICE code
clc;close all;clear;format shortEng;

Cap=1e-9;
Res=50;
Iampsrc=1;
Vampsrc=1;
Resp=1; % parallel resistence to current source
freq=5.2e9;
start_time=0; end_time=5e-9;
delt=(1/freq)/10;
Vnode2pre=0;
index=0;
Isat=5e-6; % 5uA
Vth=0.026; % 26 mV
Vdiode=0;
Vdiodepre=Vdiode;
LHS=[0;0;0;];
LHSpre=LHS;
tic;
Idiodepre=Isat*(exp(Vdiode/Vth)-1);
G0=(Isat/Vth)*exp(Vdiodepre/Vth);
I0=Idiodepre-G0*Vdiodepre;
   % MNA matrix
%     MNA=[G0+(1/Resp)         -G0             1;...
%          -G0        (1/Res)+G0      0;...
%           0          Cap/delt      -1;];  

MNA=[G0        -G0             2;...
         -G0        (1/Res)+G0      0;...
          1          Cap/delt      -1;];  

%  MNA=[G0+(1/Resp)         -G0;...
%          -G0        (1/Res)+G0+(Cap/delt);];  
for time=start_time:delt:end_time
    %Vdiodepre=0;    
    
    Isrc=0;%Iampsrc*sin(2*pi*freq*time);
    Vsrc=Vampsrc*sin(2*pi*freq*time);
 
    tol=1; %tolerance for stopping the iteration
    while(tol>1e-5)
        % RHS vector
        RHS=[Isrc-I0;...
            I0;...
            (Cap/delt)*Vnode2pre+Vsrc];

%         RHS=[Isrc-I0;...
%             I0;];
        % LHS is the solution vector
        LHS=MNA\RHS;
        Vnode2pre=LHS(2)-0;  
        Vdiode=(LHS(1)-LHS(2));   
        Vdiodepre=Vdiode;
        Idiodepre=Isat*(exp(Vdiode/Vth)-1);        
        % G0 and I0 are modified
        G0=(Isat/Vth)*exp(Vdiodepre/Vth);  
        I0=Idiodepre-G0*Vdiodepre; 
        %MNA matrix is updated
%                 MNA=[G0+(1/Resp)         -G0             0;...
%                     -G0         (1/Res)+G0      1;...
%                     0          Cap/delt         -1;];


     MNA=[G0        -G0             2;...
         -G0        (1/Res)+G0      0;...
          1          Cap/delt      -1;];  


%         MNA=[G0+(1/Resp)         -G0;...
%          -G0        (1/Res)+G0+(Cap/delt);];  
        % computing the tolerance
        tol=sum(abs([LHS(1)-LHSpre(1);LHS(2)-LHSpre(2);]));
        % storing the value of the voltage across diode
        
        LHSpre=LHS;
    end
    
    index=index+1;
    V1(index)=LHS(1);V2(index)=LHS(2);
end
toc;
figure(1),plot(0:delt:end_time,V1);grid on;figure(2),plot(0:delt:end_time,V2);grid on;
figure(3);plot(0:delt:end_time,V1-V2);grid on;
end