function MyFunctionSPICE_diode_series_RC_nw_ver6
% written on 04-May-2022 to verify the SPICE code with LTspice or ADS.
% Improved diode model is used here.
% Trying out with model given in the ADS or SPICE help file.
% Current source is used instead of a voltage source. Reverse breakdown current
% contribution is considered.
% DATE:- 11-May-2022
% AUTHOR:- Athul Shaji

clc;clear;format long;%close all;

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
Rsrc=1e20;
freq=5.2e9;
% load
Cload=1e-9;
Rload=850;
% Matching n/w
Cmatch=172.514e-15;
Lmatch=3.51485e-9;
% simulation parameters
start_time=0; end_time=5e-9;
delt=1e-12;%(1/freq)/100;
delt_old=delt;
Vnode12pre=0;
Vnode34pre=0;
Vnode40pre=0;
ILmatchpre=0;

index=0;

if(1)
    time_slice=1000;
    delt=delt/time_slice;
end

V1=zeros(length(start_time:delt_old:end_time),1);
V2=zeros(length(start_time:delt_old:end_time),1);

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
condition=1;
tic;
for time=start_time:delt_old:end_time      
    %tic
    Isrc=Iampsrc*sin(2*pi*freq*time);
    tol=1e9; %tolerance for stopping the iteration   
    for n=time:delt:time+delt_old-delt
        tol=1e9;
        Vnode12pre=LHS(1)-LHS(2);
        Vnode40pre=LHS(4)-0;
        ILmatchpre=LHS(8);        
        Vnode34pre=LHS(3)-LHS(4);
        Vdiode=(LHS(3)-LHS(4));
        Vdiodepre=Vdiode;
        while(tol>tol_value)
            % RHS vector
            RHS=[Isrc;...
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
            %tol=norm(LHS(2:4)-LHSpre(2:4)).^2;
            % storing the value of the LHS to LHSpre
            LHSpre=LHS;%Vdiodepre=Vdiode;
        end    
    end
    index=index+1;
    V0(index)=LHS(1);
    V1(index)=LHS(2);V2(index)=LHS(4);    
    diode_current(index)=Idiodepre;
    %toc
end
toc;
figure(1),plot(0:delt_old:end_time,V1);grid on;hold on;figure(2),plot(0:delt_old:end_time,V2);grid on;hold on;

Vleft=V1;Vright=V2;Idiode=diode_current;
My_MATLAB_Conversion_Efficiency(Vleft,Vright,Idiode,delt_old,freq)
end

function My_MATLAB_Conversion_Efficiency(Vleft,Vright,Idiode,delt,freq)
Vleft=Vleft(1:end);
Vright=Vright(1:end);
Idiode=Idiode(1:end);
NN=1e7;
Vleft=[Vleft,zeros(1,NN)];
Idiode=[Idiode,zeros(1,NN)];
Vright=[Vright,zeros(1,NN)];
LL=length(Vleft);
START=1;
STOP=length(Vleft);
INPUT=Vleft(START:STOP); %(V_antenna_port/dely).^2/(2*120*pi);
% Input
fft_Vin=fftshift(fft(INPUT))/LL;
%fft_Vin=fftshift(fft(V_txline_sample));
LL=length(fft_Vin);
ff=(1/delt)*(-LL/2:LL/2-1)/LL;
%figure(1),hold on;plot(ff/1e9,log10(abs(fft_Vin)));grid on;xlim([-10,10])
fft_Iin=fftshift(fft(Idiode(START:STOP)))/LL;
LL=length(fft_Iin);
ff=(1/delt)*(-LL/2:LL/2-1)/LL;
%figure(2),hold on;plot(ff/1e9,(abs(fft_V22)));grid on;xlim([-10,10])
Pin=0.5*real(fft_Vin .* conj(fft_Iin));
%Pin=abs(fft_Vin.^2)/50;
 %figure(3),hold on;plot(ff/1e9,(abs(Pin)));grid on;xlim([-30,30])

% Output
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


% eff calc
pin_val=abs(interp1(ff,(Pin),freq*0,'nearest')) +...
    2*(abs(interp1(ff,(Pin),freq*1,'nearest')) +...
    abs(interp1(ff,(Pin),freq*2,'nearest')) +...
    abs(interp1(ff,(Pin),freq*3,'nearest')) +...
    abs(interp1(ff,(Pin),freq*4,'nearest')));
%pin_val=2*(abs(interp1(ff,(Pin),2.4e9,'nearest')));
pout_val=abs(interp1(ff,(Pout),0,'nearest'));

%pin_val=(10^(-10/10))/1000;

eff=(pout_val/pin_val)*100;
disp(eff);

end