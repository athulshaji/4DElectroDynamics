function MyDiodeModel
clc;clear;close all;
index = 1;
fid = fopen('diode.txt');
while(~feof(fid))
    fgetl(fid);
    str = fgetl(fid);
    m = str2double(regexp(str, '	', 'split'));
    V_diode(index) = m(1);
    I_diode(index) = m(3);
    index = index+1;
    if(~feof(fid))
        fgetl(fid);
    end
end
fclose(fid);
figure, plot(V_diode, I_diode, 'LineWidth', 3);
axis tight, grid on;
save V_diode.mat V_diode
save I_diode.mat I_diode
% n = 1.1;
% Vth = 26e-3; %volts
% V = 0:0.001:0.8;
% I0 = 1e-12; % A/cm2
% I = I0*(exp(V./(n*Vth))-1);
% figure, plot(V,I,'LineWidth',3);
% axis tight; grid on
end