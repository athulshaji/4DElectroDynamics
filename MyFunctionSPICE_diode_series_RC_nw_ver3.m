function MyFunctionSPICE_diode_series_RC_nw_ver3
% written on 11-Apr-2022 to verify the SPICE code. Improved diode model.
% Trying out with model given in the ADS or SPICE help file
clc;clear;format shortEng;%close all;

% diode parameters
Rs=20;
Isat=5e-6; % 5uA
Vth=26e-3; % 25.85 mV @ 300K
Vdiode=0;
Vdiodepre=Vdiode;
N=1.05;
Tt=1e-11;
Cj0=0.14e-12;
Vj=0.34;
M=0.4;
Fc=0.5;
Bv=2;
Ibv=1e-4;
Xti=2;
Eg=0.69;
AREA=1;
tol_value=1e-6;
LHS=[0;0;0;0;0;0;0];
LHSpre=LHS;
if Vdiode<(-10*N*Vth)
    G0=(Isat/(N*Vth))*exp(-10);
    Idiodepre=(Isat*(exp(-10)-1)+G0*(Vdiodepre+10*N*Vth));    
else
    G0=(Isat/(N*Vth))*exp(Vdiodepre/(N*Vth));
    Idiodepre=Isat*(exp(Vdiodepre/(N*Vth))-1);
end
I0=Idiodepre-G0*Vdiodepre;
% source parameters
Iampsrc=1;
Vampsrc=1;
Rser=100;
freq=5.2e9;
% load
Cload=1e-6;
Rload=500;
% simulation parameters
start_time=0; end_time=5e-9;
delt=1e-12;%(1/freq)/100;
Vnode3pre=0;
Vnode40pre=0;
Vnode34pre=0;
index=0;

MNA=zeros(7);
MNA(1,1)=(1/Rser);
MNA(1,2)=(-1/Rser);
MNA(2,1)=(-1/Rser);
MNA(2,2)=(1/Rser)+(1/Rs);
MNA(2,3)=(-1/Rs);
MNA(3,2)=(-1/Rs);
MNA(3,3)=(1/Rs)+G0;
MNA(3,4)=-G0;
MNA(4,3)=-G0;
MNA(4,4)=G0+(1/Rload);
MNA(6,4)=Cload/delt;
MNA(6,6)=-1;
MNA(4,6)=1;
MNA(5,1)=1;
MNA(1,5)=1;

if Vdiode<=Fc*Vj
    Cj=AREA*Cj0*(1-Vdiode/Vj).^-M;
elseif Vdiode>Fc*Vj
    Cj=AREA*(Cj0/(1-Fc)^M)*(1+(M/(Vj*(1-Fc)))*(Vdiode-Fc*Vj));
end

Cdiff=Tt*G0;
Cdiode=Cdiff+Cj;
MNA(3,7)=1;
MNA(4,7)=-1;
MNA(7,3)=Cdiode/delt;
MNA(7,4)=-Cdiode/delt;
MNA(7,7)=-1;

tic;
for time=start_time:delt:end_time   
    Vsrc=Vampsrc*sin(2*pi*freq*time);
    tol=1e9; %tolerance for stopping the iteration
    while(tol>tol_value)
        % RHS vector
        RHS=[0;...
            0;...
            -I0;...
            I0;...
            Vsrc;...
            (Cload/delt)*Vnode40pre;...
            (Cdiode/delt)*Vnode34pre];
        % LHS is the solution vector
        LHS=MNA\RHS;
        Vnode40pre=LHS(4)-0;  
        Vnode34pre=LHS(3)-LHS(4);
        Vdiode=(LHS(3)-LHS(4));   
        Vdiodepre=Vdiode;
         % G0 and I0 are modified
        if Vdiodepre<(-10*N*Vth)
            G0=(Isat/(N*Vth))*exp(-10);
            Idiodepre=(Isat*(exp(-10)-1)+G0*(Vdiodepre+10*N*Vth));
        else
            G0=(Isat/(N*Vth))*exp(Vdiodepre/(N*Vth));
            Idiodepre=Isat*(exp(Vdiodepre/(N*Vth))-1);
        end          
        I0=Idiodepre-G0*(LHS(3)-LHS(4));         
        % junction capacitance is evaluated             
        if Vdiodepre<=Fc*Vj
            Cj=AREA*Cj0*(1-(Vdiodepre)/Vj).^-M;
        elseif Vdiodepre>Fc*Vj
            Cj=AREA*(Cj0/(1-Fc)^M)*(1+(M/(Vj*(1-Fc)))*((Vdiodepre)-Fc*Vj));
        end
        
        Cdiff=Tt*G0;
        Cdiode=Cdiff+Cj;        
        % MNA matrix is updated
        MNA(3,3)=(1/Rs)+G0;
        MNA(3,4)=-G0;
        MNA(4,3)=-G0;
        MNA(4,4)=G0+(1/Rload);        
        MNA(7,3)=Cdiode/delt;
        MNA(7,4)=-Cdiode/delt;        
        % computing the tolerance
        %tol=sum(abs([LHS(2)-LHSpre(2);LHS(3)-LHSpre(3);LHS(4)-LHSpre(4);]).^2);
        tol=sum([(LHS(2:4)-LHSpre(2:4));].^2);
        % storing the value of the LHS to LHSpre       
        LHSpre=LHS;
    end    
    index=index+1;
    diode_capacitance(index)=Cdiode;
    V1(index)=LHS(2);V2(index)=LHS(4);
end
toc;
figure(1),plot(0:delt:end_time,V1);grid on;hold on;figure(2),plot(0:delt:end_time,V2);grid on;hold on;
%figure(3);plot(0:delt:end_time,V1-V2);grid on;
end