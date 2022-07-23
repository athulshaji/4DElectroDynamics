function MyFunctionSPICE_diode_series_RC_nw_ver7
% written on 04-May-2022 to verify the SPICE code with LTspice or ADS.
% Improved diode model is used here.
% Trying out with model given in the ADS or SPICE help file.
% Current source is used instead of a voltage source. Reverse breakdown current
% contribution is considered. Rs (diode resistance is included in the
% equation
% DATE:- 11-May-2022
% AUTHOR:- Athul Shaji

clc;clear;format shortEng;close all;

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
AREA=1;%default value
ExplI=1;%default value
Ib0=Ibv*exp(-Bv/(N*Vth));
tol_value=1e-9;
LHS=[0;0;0;0;0;0;0;0];
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
Iampsrc=1e-3;
Vampsrc=1;
Rshunt=1;
freq=5.2e9;
% load
Cload=1e-9;
Rload=50;
% Matching n/w
Cmatch=172.514e-15;
Lmatch=3.51485e-9;
% simulation parameters
start_time=0; end_time=5e-9;
delt=1e-12;%(1/freq)/100;
delt_old=delt;
Vnode12pre=0;
Vnode23pre=0;
Vnode30pre=0;
ILmatchpre=0;

index=0;

if(1)
    time_slice=1;
    delt=delt/time_slice;
end

V1=zeros(length(start_time:delt_old:end_time),1);
V2=zeros(length(start_time:delt_old:end_time),1);

MNA=zeros(7);
MNA(2,2)=G0;
MNA(2,3)=-G0;
MNA(3,2)=-G0;
MNA(3,3)=G0+(1/Rload);
MNA(4,1)=Cmatch/delt;
MNA(4,2)=-Cmatch/delt;
MNA(4,4)=-1;
MNA(2,4)=-1;
MNA(1,4)=1;

if Vdiode<=Fc*Vj
    Cj=AREA*Cj0*(1-Vdiode/Vj).^-M;
elseif Vdiode>Fc*Vj
    Cj=AREA*(Cj0/(1-Fc)^M)*(1+(M/(Vj*(1-Fc)))*(Vdiode-Fc*Vj));
end

Cdiff=Tt*G0;
Cdiode=Cdiff+Cj;
MNA(2,5)=1;
MNA(3,5)=-1;
MNA(5,2)=Cdiode/delt;
MNA(5,3)=-Cdiode/delt;
MNA(5,5)=-1;

MNA(6,3)=Cload/delt;
MNA(6,6)=-1;
MNA(3,6)=1;

MNA(7,2)=1;
MNA(7,7)=-Lmatch/delt;
MNA(2,7)=1;
func=@(Idiode) Isat*(exp((Vdiodepre-Rs*Idiode)/(N*Vth))-1)-Idiode;
options = optimoptions('fsolve','Display','off');
tic;
for time=start_time:delt_old:end_time
    %     Vsrc=Vampsrc*sin(2*pi*freq*time);
    Isrc=Iampsrc*sin(2*pi*freq*time);
    tol=1e9; %tolerance for stopping the iteration   
    for n=time:delt:time+delt_old-delt
        tol=1e9;
        while(tol>tol_value)             
            % RHS vector
            RHS=[Isrc;...                   
                -I0;...
                I0;...
                (Cmatch/delt)*Vnode12pre;...
                (Cdiode/delt)*Vnode23pre;...
                (Cload/delt)*Vnode30pre;...
                (-Lmatch/delt)*ILmatchpre;];
            % LHS is the solution vector
            LHS=MNA\RHS;
            Vnode12pre=LHS(1)-LHS(2);
            Vnode23pre=LHS(2)-LHS(3);
            Vnode30pre=LHS(3)-0;
            ILmatchpre=LHS(7);
            Vdiode=(LHS(2)-0);
            Vdiodepre=Vdiode;
            Vbmax=N*Vth*log(ExplI/Ibv);
            Gbmax=ExplI/(N*Vth);
            % G0 and I0 are modified
            if Vdiodepre<(-10*N*Vth)
                G0=(Isat/(N*Vth))*exp(-10);
                Idiodepre=(Isat*(exp(-10)-1)+G0*(Vdiodepre+10*N*Vth));
            elseif -(Vdiodepre+Bv)>Vbmax
                Idiodepre=-(ExplI+(-(Vdiodepre+Bv)-Vbmax)*Gbmax-Ib0);
                G0=Gbmax;
            else
                G0=(Isat/(N*Vth))*(exp((Vdiodepre-Rs*Idiodepre)/(N*Vth)));
                Idiodepre=fsolve(func,Idiodepre,options);                
            end
            I0=Idiodepre-G0*(LHS(2)-0);
            % junction capacitance is evaluated
            if Vdiodepre<=Fc*Vj
                Cj=AREA*Cj0*(1-(Vdiodepre)/Vj).^-M;
            elseif Vdiodepre>Fc*Vj
                Cj=AREA*(Cj0/(1-Fc)^M)*(1+(M/(Vj*(1-Fc)))*((Vdiodepre)-Fc*Vj));
            end

            Cdiff=Tt*G0;
            Cdiode=Cdiff+Cj;
            % MNA matrix is updated
            MNA(2,2)=G0;
            MNA(2,3)=-G0;
            MNA(3,2)=-G0;
            MNA(3,3)=G0+(1/Rload);
            MNA(5,2)=Cdiode/delt;
            MNA(5,3)=-Cdiode/delt;
            % computing the tolerance            
            %tol=abs(LHS(2:3)-LHSpre(2:3)).^4;
            tol=abs((LHS(2)-LHS(3))-(LHSpre(2)-LHSpre(3))).^2;
            % storing the value of the LHS to LHSpre
            LHSpre=LHS;
        end    
    end
    index=index+1;
    diode_capacitance(index)=Cdiode;
    V1(index)=LHS(2);V2(index)=LHS(3);
end
toc;
figure(1),plot(0:delt_old:end_time,V1);grid on;hold on;figure(2),plot(0:delt_old:end_time,V2);grid on;hold on;
%figure(3);plot(0:delt:end_time,V1-V2);grid on;
end