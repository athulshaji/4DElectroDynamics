function TxLine_CNFDTD_ver6_11May2022

%% Introduction
% This code combines one dimensional Crank-Nicolson Finite-Difference 
% Time-Domain (CN-FDTD-1D) with SPICE. The transmission line is
% modelled using the transmission line equations. The per unit length
% parameters are manually defined.

%% clearing the workspace, defining the format for the numbers display, closing the figures , and clearing the command window
clc;clear;format long;close all;

%% Generate the movie. The video is NOT saved.
MOVIE=0;

%% CN-FDTD-SPICE simulation setup

S=5; % stability factor. As we are using the CN-FDTD algorithm it is possible for us to go beyond the time limit given by CFl condition.
Nz=5001; % number of spatial divisions
R=0;L=2500e-9;G=0;C=1e-9; % per unit length parametes for a Lossless tx line
Lz=1.44e-2; % length of the tx line
Tsteps=5000;
freq=5.2e9; % frequency of the sinusoidal signal
AMP=0.018; 
t=0:Tsteps-1;
delz=Lz/(Nz-1);
vp=1/sqrt(L*C);
delt_cfl=delz/vp; % CFL criterion
delt=delt_cfl*S;
delt_old=delt;
simulation_time=Tsteps*delt;
VSRC(t+1)=AMP*sin(2*pi*freq*t*delt);
disp(['Simulation time: ',num2str(simulation_time),' seconds.']);


%% CN-FDTD (LHS and RHS) matrix setup
Mz=Nz-2;
Dzb=(1/(2*delz))*spdiags([-ones(Mz,1),ones(Mz,1)],[-1,0],Mz,Mz);
Dzf=(1/(2*delz))*spdiags([-ones(Mz,1),ones(Mz,1)],[0,1],Mz,Mz);

LHS_1=(C/delt)*speye(Mz);
LHS_2=Dzb;
LHS_3=Dzf;
LHS_4=(L/delt)*speye(Mz);

RHS_1=(C/delt)*speye(Mz);
RHS_2=-Dzb;
RHS_3=-Dzf;
RHS_4=(L/delt)*speye(Mz);

LHS=[LHS_1 LHS_2;LHS_3 LHS_4];
tic;[L_,U_,P_]=lu(LHS);toc;

RHS=[RHS_1 RHS_2;RHS_3 RHS_4];

RHS_additional=[zeros(Mz,1);zeros(Mz,1)];
% the voltage and current as a SINGLE vector: vector of unknowns
VI=[zeros(Mz,1);zeros(Mz,1)];

%% The diode parameters used for the SPICE modelling is provided
% below. These values are provided in the data sheet
Rs=20;
Isat=5e-6; % 5uA
Vth=25.85e-3; % 25.85 mV @ 300K
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
ExplI=2;%default value
Ib0=Ibv*exp(-Bv/(N*Vth));
tol_value=1e-3;
LHSckt=[0;0;0;0;0;0;0;0];
LHScktpre=LHSckt;
if Vdiode<(-10*N*Vth)
    G0=(Isat/(N*Vth))*exp(-10);
    Idiodepre=(Isat*(exp(-10)-1)+G0*(Vdiodepre+10*N*Vth));
else
    G0=(Isat/(N*Vth))*exp(Vdiodepre/(N*Vth));
    Idiodepre=Isat*(exp(Vdiodepre/(N*Vth))-1);
end
I0=Idiodepre-G0*Vdiodepre;
% load
Cload=1e-12;
Rload=500;
% Matching n/w
Cmatch=172.514e-15;
Lmatch=3.51485e-9;

% Here this time slice factor is used to discretize the time step provided
% so that the time resolution can be increased. So if the original time
% step for the source is 1ns then the simulation within each of the time
% steps is evaluated at 1ps
if(1)
    time_slice=1000;
    delt=delt/time_slice;
end

%% Initializing the SPICE
[MNAckt,LHSckt,LHScktpre,G0,I0,Vdiode,Cdiode,Cmatch,Lmatch,Rload,Cload]=SPICE_Init(delt);RHSckt=LHSckt;

% Intializing the variables
Ickt=0;index=0;V0=[];V1=[];V2=[];I_branch=0;V2pre=0;
Vnode12pre=0;
Vnode34pre=0;
Vnode40pre=0;
ILmatchpre=0;
input_source=0;

%% CNFDTD + SPICE section
tic
for n=1:Tsteps
    RHS_additional(1)=S*(C/delt_old)*(-VSRC(n));
    VI=(U_\(L_\(P_*(RHS*VI+RHS_additional))));
    input_source=-(VI(2*Mz-1)+Ickt)/2;
    % SPICE stitching with the CN-FDTD
    time=(n-1)*delt_old;
    [MNAckt,RHSckt,LHSckt,LHScktpre,V0,V1,V2,G0,I0,Vnode12pre,Vnode40pre,ILmatchpre,Vdiode,Vdiodepre,tol_value,index,I_branch,diode_current,src_current]=MySPICE(MNAckt,RHSckt,LHSckt,LHScktpre,Cdiode,Cload,Rload,Cmatch,Lmatch,V0,V1,V2,G0,I0,Vnode12pre,Vnode40pre,ILmatchpre,Vdiode,Vdiodepre,tol_value,index,I_branch,input_source,time,delt,delt_old);
    % updating the current value in the last node of the tx-line using the
    % voltage generated across the
    % SPICE circuit
    RHS_additional(2*Mz)=-((V1(n)+V2pre))/(2*delz);
    Ickt=I_branch;
    V2pre=V1(n);

    if (MOVIE)
        plot(VI(1:Mz));grid on;title(num2str(n));ylim([-2.1,2.1]);drawnow;
    end

% Saving the data at the end of the simulation
    if(mod(n,Tsteps)==0)
        save case-S=5.mat
    end


end
toc

Vleft=V1;Vright=V2;Idiode=diode_current;Isrc=src_current;
% My MATLAB function for evaluating the RF-DC conversion efficiency.
% Function definition can be seen below.
My_MATLAB_Conversion_Efficiency(Vleft,Isrc,Vright,Idiode,delt_old,freq)
end


function [MNA,LHS,LHSpre,G0,I0,Vdiode,Cdiode,Cmatch,Lmatch,Rload,Cload]=SPICE_Init(delt)
% diode parameters
Rs=20;
Isat=5e-6; % 5uA
Vth=25.85e-3; % 25.85 mV @ 300K
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
ExplI=2;%default value
Ib0=Ibv*exp(-Bv/(N*Vth));
tol_value=1e-3;
LHS=[0;0;0;0;0;0;0;0;];
LHSpre=LHS;
if Vdiode<(-10*N*Vth)
    G0=(Isat/(N*Vth))*exp(-10);
    Idiodepre=(Isat*(exp(-10)-1)+G0*(Vdiodepre+10*N*Vth));
else
    G0=(Isat/(N*Vth))*exp(Vdiodepre/(N*Vth));
    Idiodepre=Isat*(exp(Vdiodepre/(N*Vth))-1);
end
I0=Idiodepre-G0*Vdiodepre;
% load
Cload=1e-12;
Rload=500;
% Matching n/w
Cmatch=172.514e-15;
Lmatch=3.51485e-9;
Csrc=1e-9;
MNA=zeros(8);
%MNA(1,1)=1/Rsrc;
MNA(2,2)=1/Rs;
MNA(2,3)=-1/Rs;
MNA(3,2)=-1/Rs;
MNA(3,3)=(1/Rs)+G0;
MNA(3,4)=-G0;
MNA(4,3)=-G0;
MNA(4,4)=G0+(1/Rload);
MNA(5,1)=Cmatch/delt;
MNA(5,2)=-Cmatch/delt;
MNA(5,5)=-1;
MNA(2,5)=-1;
MNA(1,5)=1;
%Capacitor
% MNA(9,1)=1;
% MNA(1,9)=1;
% MNA(9,9)=-1;


if Vdiode<=Fc*Vj
    Cj=AREA*Cj0*(1-Vdiode/Vj).^-M;
elseif Vdiode>Fc*Vj
    Cj=AREA*(Cj0/(1-Fc)^M)*(1+(M/(Vj*(1-Fc)))*(Vdiode-Fc*Vj));
end

Cdiff=Tt*G0;
Cdiode=Cdiff+Cj;
MNA(3,6)=1;
MNA(4,6)=-1;
MNA(6,3)=Cdiode/delt;
MNA(6,4)=-Cdiode/delt;
MNA(6,6)=-1;

MNA(7,4)=Cload/delt;
MNA(7,7)=-1;
MNA(4,7)=1;

MNA(8,2)=1;
MNA(8,8)=-Lmatch/delt;
MNA(2,8)=1;
end

function [MNA,RHS,LHS,LHSpre,V0,V1,V2,G0,I0,Vnode12pre,Vnode40pre,ILmatchpre,Vdiode,Vdiodepre,tol_value,index,I_branch,diode_current,src_current]=MySPICE(MNA,RHS,LHS,LHSpre,Cdiode,Cload,Rload,Cmatch,Lmatch,V0,V1,V2,G0,I0,Vnode12pre,Vnode40pre,ILmatchpre,Vdiode,Vdiodepre,tol_value,index,I_branch,input_source,time,delt,delt_old,tolpre)
Rs=20;
Isat=5e-6; % 5uA
Vth=25.85e-3; % 25.85 mV @ 300K
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
ExplI=2;%default value
Ib0=Ibv*exp(-Bv/(N*Vth));
Csrc=1e-9;
for n=time:delt:time+delt_old-delt
    tol=1e9;
    SRC=input_source;
    Vnode10pre=LHS(1)-0;
    Vnode12pre=LHS(1)-LHS(2);
    Vnode40pre=LHS(4)-0;
    ILmatchpre=LHS(8);
    Vnode34pre=LHS(3)-LHS(4);
    Vdiode=(LHS(3)-LHS(4));
    Vdiodepre=Vdiode;
    tolpre=1e9;
    while(tol>tol_value)
        % RHS vector
        RHS=[SRC;...
            0;...
            -I0;...
            I0;...
            (Cmatch/delt)*Vnode12pre;...
            (Cdiode/delt)*Vnode34pre;...
            (Cload/delt)*Vnode40pre;...
            (-Lmatch/delt)*ILmatchpre;];

        % LHS is the solution vector
        LHS=MNA\RHS;
        Vdiode=(LHS(3)-LHS(4));%Vnode34pre=LHS(3)-LHS(4);
        Vdiodepre=Vdiode;
        Vmax=N*Vth*log((ExplI/Isat)+1);
        Vbmax=N*Vth*log(ExplI/Ibv);
        Gbmax=ExplI/(N*Vth);
        % G0 and I0 are modified using applicable equations at different voltage ranges
        if Vdiodepre<(-10*N*Vth)
            if -(Bv+Vbmax)>Vdiodepre
                Idiodepre=-(ExplI+(-((Vdiodepre)+Bv)-Vbmax)*Gbmax-Ib0);
                G0=Gbmax;
            elseif -Bv>Vdiodepre
                Idiodepre=-Ibv*exp(-((Vdiodepre)+Bv)/(N*Vth))+Ib0;
                G0=-Idiodepre/(N*Vth);
            else
                G0=(Isat/(N*Vth))*exp(-10);
                Idiodepre=(Isat*(exp(-10)-1)+G0*(Vdiodepre+10*N*Vth));
            end
        elseif (-10*N*Vth<=Vdiodepre && Vdiodepre<=Vmax)
            G0=(Isat/(N*Vth))*exp(Vdiodepre/(N*Vth));
            Idiodepre=Isat*(exp(Vdiodepre/(N*Vth))-1);
        else
            Idiodepre=0;
            G0=0;
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
        MNA(6,3)=Cdiode/delt;
        MNA(6,4)=-Cdiode/delt;
        % computing the tolerance
        tol=abs((LHS(3)-LHS(4))-(LHSpre(3)-LHSpre(4)));
%             if(abs(tol-tolpre)<1e-1)
%                 break;
%             else 
%                 tolpre=tol;
%             end
        %tol=norm(LHS(2:4)-LHSpre(2:4)).^2;
        % storing the value of the LHS to LHSpre
        LHSpre=LHS;%Vdiodepre=Vdiode;
    end
end
index=index+1;
V0(index)=LHS(1);
V1(index)=LHS(2);V2(index)=LHS(4);
I_branch=SRC;
diode_current(index)=Idiodepre;
src_current(index)=SRC;
    %toc
end

function My_MATLAB_Conversion_Efficiency(Vleft,Isrc,Vright,Idiode,delt,freq)
Vleft=Vleft(1:end);
Vright=Vright(1:end);
Idiode=Idiode(1:end);
Isrc=Isrc(1:end);
NN=1e6; % appending the data with zeros for finer freq resolution
Vleft=[Vleft,zeros(1,NN)];
Idiode=[Idiode,zeros(1,NN)];
Isrc=[Isrc,zeros(1,NN)];
Vright=[Vright,zeros(1,NN)];
LL=length(Vleft);
START=2500;
STOP=length(Vleft);
INPUT=Vleft(START:STOP); %(V_antenna_port/dely).^2/(2*120*pi);

%% Input
fft_Vin=fftshift(fft(INPUT))/LL;
%fft_Vin=fftshift(fft(V_txline_sample));
LL=length(fft_Vin);
ff=(1/delt)*(-LL/2:LL/2-1)/LL;
%figure(1),hold on;plot(ff/1e9,log10(abs(fft_Vin)));grid on;xlim([-10,10])
fft_Iin=fftshift(fft(Isrc(START:STOP)))/LL;
LL=length(fft_Iin);
ff=(1/delt)*(-LL/2:LL/2-1)/LL;
%figure(2),hold on;plot(ff/1e9,(abs(fft_V22)));grid on;xlim([-10,10])
Pin=0.5*real(fft_Vin .* conj(fft_Iin));
%Pin=abs(fft_Vin.^2)/50;
 %figure(3),hold on;plot(ff/1e9,(abs(Pin)));grid on;xlim([-30,30])

%% Output
fft_Vout=fftshift(fft(Vright(START:STOP)))/LL;
LL=length(fft_Vout);
ff=(1/delt)*(-LL/2:LL/2-1)/LL;
%figure(4),hold on;plot(ff/1e9,log10(abs(fft_Vout)));grid on;xlim([-10,10])
fft_Iout=fftshift(fft(Idiode(START:STOP)))/LL;
LL=length(fft_Iout);
ff=(1/delt)*(-LL/2:LL/2-1)/LL;
%figure(5),hold on;plot(ff/1e9,(abs(fft_V22)));grid on;xlim([-10,10])
Pout=0.5*real(fft_Vout .* conj(fft_Iout));
 %figure(6),hold on;plot(ff/1e9,(abs(Pout)));grid on;xlim([-30,30])
%figure(10),hold on;plot(V_antenna_port);


% calculating the input power (RF power)
pin_val=abs(interp1(ff,(Pin),freq*0,'nearest')) +...
    2*(abs(interp1(ff,(Pin),freq*1,'nearest')) +...
    abs(interp1(ff,(Pin),freq*2,'nearest')) +...
    abs(interp1(ff,(Pin),freq*3,'nearest')) +...
    abs(interp1(ff,(Pin),freq*4,'nearest')));

% calculating the output power (DC power)
pout_val=abs(interp1(ff,(Pout),0,'nearest'));

% efficiency calculation
eff=(pout_val/pin_val)*100;
disp(['RF-DC conversion efficiency: ',num2str(eff),'%']);
end

