clc;close all;clear;
% without diode
load('D:\OneDrive - Indian Institute of Science\Dev\MATLAB_Code\MATLAB_Code-master\new_codes\02Oct2021\case10.mat')
save('02Oct2021\data\delt.txt','delt','-ascii');
save('02Oct2021\data\VC1_1.txt','VC1','-ascii');
% near-field
save('02Oct2021\data\nffieldxx_1.txt','nffieldxx','-ascii');
save('02Oct2021\data\nffieldyy_1.txt','nffieldyy','-ascii');
save('02Oct2021\data\nffieldzz_1.txt','nffieldzz','-ascii');
% far-field
save('02Oct2021\data\fffieldxx_1.txt','fffieldxx','-ascii');
save('02Oct2021\data\fffieldyy_1.txt','fffieldyy','-ascii');
save('02Oct2021\data\fffieldzz_1.txt','fffieldzz','-ascii');

clear
% with diode
load('D:\OneDrive - Indian Institute of Science\Dev\MATLAB_Code\MATLAB_Code-master\new_codes\02Oct2021\case11.mat')
save('02Oct2021\data\VC1_2.txt','VC1','-ascii');
% near-field
save('02Oct2021\data\nffieldxx_2.txt','nffieldxx','-ascii');
save('02Oct2021\data\nffieldyy_2.txt','nffieldyy','-ascii');
save('02Oct2021\data\nffieldzz_2.txt','nffieldzz','-ascii');
% far-field
save('02Oct2021\data\fffieldxx_2.txt','fffieldxx','-ascii');
save('02Oct2021\data\fffieldyy_2.txt','fffieldyy','-ascii');
save('02Oct2021\data\fffieldzz_2.txt','fffieldzz','-ascii');

clear
% case 1: without diode
load('delt.txt');
load('02Oct2021\data\VC1_1.txt');
load('02Oct2021\data\nffieldxx_1.txt');
load('02Oct2021\data\nffieldyy_1.txt');
load('02Oct2021\data\nffieldzz_1.txt');
load('02Oct2021\data\fffieldxx_1.txt');
load('02Oct2021\data\fffieldyy_1.txt');
load('02Oct2021\data\fffieldzz_1.txt');
% case 2: with diode
load('02Oct2021\data\VC1_2.txt');
load('02Oct2021\data\nffieldxx_2.txt');
load('02Oct2021\data\nffieldyy_2.txt');
load('02Oct2021\data\nffieldzz_2.txt');
load('02Oct2021\data\fffieldxx_2.txt');
load('02Oct2021\data\fffieldyy_2.txt');
load('02Oct2021\data\fffieldzz_2.txt');

% padding the zeros after 3500 samples
VC1_1(3501:end)=0;VC1_2(3501:end)=0;

nffieldxx_1(3501:end)=0;nffieldxx_2(3501:end)=0;
nffieldyy_1(3501:end)=0;nffieldyy_2(3501:end)=0;
nffieldzz_1(3501:end)=0;nffieldzz_2(3501:end)=0;

fffieldxx_1(3501:end)=0;fffieldxx_2(3501:end)=0;
fffieldyy_1(3501:end)=0;fffieldyy_2(3501:end)=0;
fffieldzz_1(3501:end)=0;fffieldzz_2(3501:end)=0;

% frequency range
L=length(nffieldxx_1);
f=(1/delt)*(-L/2:L/2-1)/L;

fft_VC1_1=fftshift(fft(VC1_1))/L;fft_VC1_2=fftshift(fft(VC1_2))/L;

fft_nffieldxx_1=fftshift(fft(nffieldxx_1))/L;fft_nffieldxx_2=fftshift(fft(nffieldxx_2))/L;
fft_nffieldyy_1=fftshift(fft(nffieldyy_1))/L;fft_nffieldyy_2=fftshift(fft(nffieldyy_2))/L;
fft_nffieldzz_1=fftshift(fft(nffieldzz_1))/L;fft_nffieldzz_2=fftshift(fft(nffieldzz_2))/L;

fft_fffieldxx_1=fftshift(fft(fffieldxx_1))/L;fft_fffieldxx_2=fftshift(fft(fffieldxx_2))/L;
fft_fffieldyy_1=fftshift(fft(fffieldyy_1))/L;fft_fffieldyy_2=fftshift(fft(fffieldyy_2))/L;
fft_fffieldzz_1=fftshift(fft(fffieldzz_1))/L;fft_fffieldzz_2=fftshift(fft(fffieldzz_2))/L;

% generating figures
figure(1),plot(log10(abs(VC1_1)));hold on;grid on;plot(log10(abs(VC1_2)));hold on;grid on;xlabel('time steps');ylabel('log10(|antenna terminal voltages|)');legend('antenna with load resistor','antenna-matching nw-shunt diode-tx line-load resistor');
figure(2),plot(log10(abs(VC1_2-VC1_1)));hold on;grid on;xlabel('time steps');ylabel('log10(|difference voltage at antenna terminals|)');
figure(3),plot(f,log10(abs(fft_VC1_1)));hold on;grid on;plot(f,log10(abs(fft_VC1_2)));hold on;grid on;xlim([-1e10,1e10]);xlabel('freq(Hz)');ylabel('log10(|FFT(voltage at antenna terminals)|)');legend('antenna with load resistor','antenna-matching nw-shunt diode-tx line-load resistor');
figure(4),plot(f,log10(abs(fft_VC1_2-fft_VC1_1)));hold on;grid on;xlim([-1e10,1e10]);xlim([-1e10,1e10]);xlabel('freq(Hz)');ylabel('log10(|FFT(difference of voltages at antenna terminals)|)');

figure(5),plot(log10(abs(nffieldxx_1)));hold on;grid on;plot(log10(abs(nffieldxx_2)));grid on;hold on;xlabel('time steps');ylabel('log10(|near-field Ex|)');legend('antenna with load resistor','antenna-matching nw-shunt diode-tx line-load resistor');
figure(6),plot(log10(abs(nffieldxx_2-nffieldxx_1)));hold on;grid on;xlabel('time steps');ylabel('log10(|difference of near field Ex|)');
figure(7),plot(f,log10(abs(fft_nffieldxx_1)));hold on;grid on;plot(f,log10(abs(fft_nffieldxx_2)));hold on;grid on;xlim([-1e10,1e10]);xlabel('freq(Hz)');ylabel('log10(|FFT(near-field Ex)|)');legend('antenna with load resistor','antenna-matching nw-shunt diode-tx line-load resistor');
figure(8),plot(f,log10(abs(fft_nffieldxx_2-fft_nffieldxx_1)));hold on;grid on;xlim([-1e10,1e10]);xlabel('freq(Hz)');ylabel('log10(|FFT(difference of near field Ex)|)');

figure(9),plot(log10(abs(nffieldyy_1)));hold on;grid on;plot(log10(abs(nffieldyy_2)));grid on;hold on;xlabel('time steps');ylabel('log10(|near-field Ey|)');legend('antenna with load resistor','antenna-matching nw-shunt diode-tx line-load resistor');
figure(10),plot(log10(abs(nffieldyy_2-nffieldyy_1)));hold on;grid on;xlabel('time steps');ylabel('log10(|difference of near field Ey|)');
figure(11),plot(f,log10(abs(fft_nffieldyy_1)));hold on;grid on;plot(f,log10(abs(fft_nffieldyy_2)));hold on;grid on;xlim([-1e10,1e10]);xlim([-1e10,1e10]);xlabel('freq(Hz)');ylabel('log10(|FFT(near-field Ey)|)');legend('antenna with load resistor','antenna-matching nw-shunt diode-tx line-load resistor');
figure(12),plot(f,log10(abs(fft_nffieldyy_2-fft_nffieldyy_1)));hold on;grid on;xlim([-1e10,1e10]);xlabel('freq(Hz)');ylabel('log10(|FFT(difference of near field Ey)|)');

figure(13),plot(log10(abs(nffieldzz_1)));hold on;grid on;plot(log10(abs(nffieldzz_2)));grid on;hold on;xlabel('time steps');ylabel('log10(|near-field Ez|)');legend('antenna with load resistor','antenna-matching nw-shunt diode-tx line-load resistor');
figure(14),plot(log10(abs(nffieldzz_2-nffieldzz_1)));hold on;grid on;xlabel('time steps');ylabel('log10(|difference of near field Ez|)');
figure(15),plot(f,log10(abs(fft_nffieldzz_1)));hold on;grid on;plot(f,log10(abs(fft_nffieldzz_2)));hold on;grid on;xlim([-1e10,1e10]);xlabel('freq(Hz)');ylabel('log10(|FFT(near-field Ez)|)');legend('antenna with load resistor','antenna-matching nw-shunt diode-tx line-load resistor');
figure(16),plot(f,log10(abs(fft_nffieldzz_2-fft_nffieldzz_1)));hold on;grid on;xlim([-1e10,1e10]);xlabel('freq(Hz)');ylabel('log10(|FFT(difference of near field Ez)|)');

figure(17),plot(log10(abs(fffieldxx_1)));hold on;grid on;plot(log10(abs(fffieldxx_2)));grid on;hold on;xlabel('time steps');ylabel('log10(|far-field Ex|)');legend('antenna with load resistor','antenna-matching nw-shunt diode-tx line-load resistor');
figure(18),plot(log10(abs(fffieldxx_2-fffieldxx_1)));hold on;grid on;xlabel('time steps');ylabel('log10(|difference of far field Ex|)');
figure(19),plot(f,log10(abs(fft_fffieldxx_1)));hold on;grid on;plot(f,log10(abs(fft_fffieldxx_2)));hold on;grid on;xlim([-1e10,1e10]);xlabel('freq(Hz)');ylabel('log10(|FFT(far-field Ex)|)');legend('antenna with load resistor','antenna-matching nw-shunt diode-tx line-load resistor');
figure(20),plot(f,log10(abs(fft_fffieldxx_2-fft_fffieldxx_1)));hold on;grid on;xlim([-1e10,1e10]);xlabel('freq(Hz)');ylabel('log10(|FFT(difference of far field Ex)|)');

figure(21),plot(log10(abs(fffieldyy_1)));hold on;grid on;plot(log10(abs(fffieldyy_2)));grid on;hold on;xlabel('time steps');ylabel('log10(|far-field Ey|)');legend('antenna with load resistor','antenna-matching nw-shunt diode-tx line-load resistor');
figure(22),plot(log10(abs(fffieldyy_2-fffieldyy_1)));hold on;grid on;xlabel('time steps');ylabel('log10(|difference of far field Ey|)');
figure(23),plot(f,log10(abs(fft_fffieldyy_1)));hold on;grid on;plot(f,log10(abs(fft_fffieldyy_2)));hold on;grid on;xlim([-1e10,1e10]);xlabel('freq(Hz)');ylabel('log10(|FFT(far-field Ey)|)');legend('antenna with load resistor','antenna-matching nw-shunt diode-tx line-load resistor');
figure(24),plot(f,log10(abs(fft_fffieldyy_2-fft_fffieldyy_1)));hold on;grid on;xlim([-1e10,1e10]);xlabel('freq(Hz)');ylabel('log10(|FFT(difference of far field Ey)|)');

figure(25),plot(log10(abs(fffieldzz_1)));hold on;grid on;plot(log10(abs(fffieldzz_2)));grid on;hold on;xlabel('time steps');ylabel('log10(|far-field Ez|)');legend('antenna with load resistor','antenna-matching nw-shunt diode-tx line-load resistor');
figure(26),plot(log10(abs(fffieldzz_2-fffieldzz_1)));hold on;grid on;xlabel('time steps');ylabel('log10(|difference of far field Ez|)');
figure(27),plot(f,log10(abs(fft_fffieldzz_1)));hold on;grid on;plot(f,log10(abs(fft_fffieldzz_2)));hold on;grid on;xlim([-1e10,1e10]);xlabel('freq(Hz)');ylabel('log10(|FFT(far-field Ez)|)');legend('antenna with load resistor','antenna-matching nw-shunt diode-tx line-load resistor');
figure(28),plot(f,log10(abs(fft_fffieldzz_2-fft_fffieldzz_1)));hold on;grid on;xlim([-1e10,1e10]);xlabel('freq(Hz)');ylabel('log10(|FFT(difference of far field Ez)|)');

clear;

% clc;close all;clear;
% load('D:\OneDrive - Indian Institute of Science\Dev\MATLAB_Code\MATLAB_Code-master\new_codes\02Oct2021\case8.mat')
% LL=2300;f=(1/delt)*(-LL/2:LL/2-1)/LL;
% figure(1);plot(log10(abs(fieldxx(1:LL))));grid on;hold on;
% figure(2);plot(log10(abs(fieldyy(1:LL))));grid on;hold on;
% figure(3);plot(log10(abs(fieldzz(1:LL))));grid on;hold on;
% figure(4),plot(f,log10(abs(fftshift(fft(fieldxx(1:LL)))/LL)));grid on;hold on;xlim([-1e10,1e10]);
% figure(5),plot(f,log10(abs(fftshift(fft(fieldyy(1:LL)))/LL)));grid on;hold on;xlim([-1e10,1e10]);
% figure(6),plot(f,log10(abs(fftshift(fft(fieldzz(1:LL)))/LL)));grid on;hold on;xlim([-1e10,1e10]);
% figure(7),plot(f,log10(abs(fftshift(fft(VC1(1:LL)))/LL)));grid on;hold on;xlim([-1e10,1e10]);
% figure(8),plot(f,log10(abs(fftshift(fft(VC2(1:LL)))/LL)));grid on;hold on;xlim([-1e10,1e10]);
% figure(9),plot(VC1(1:LL));grid on;hold on;
% clear
% load('D:\OneDrive - Indian Institute of Science\Dev\MATLAB_Code\MATLAB_Code-master\new_codes\02Oct2021\case9.mat')
% LL=2300;f=(1/delt)*(-LL/2:LL/2-1)/LL;
% figure(1);plot(log10(abs(fieldxx(1:LL))));grid on;hold on;
% figure(2);plot(log10(abs(fieldyy(1:LL))));grid on;hold on;
% figure(3);plot(log10(abs(fieldzz(1:LL))));grid on;hold on;
% figure(4),plot(f,log10(abs(fftshift(fft(fieldxx(1:LL)))/LL)));grid on;hold on;xlim([-1e10,1e10]);
% figure(5),plot(f,log10(abs(fftshift(fft(fieldyy(1:LL)))/LL)));grid on;hold on;xlim([-1e10,1e10]);
% figure(6),plot(f,log10(abs(fftshift(fft(fieldzz(1:LL)))/LL)));grid on;hold on;xlim([-1e10,1e10]);
% figure(7),plot(f,log10(abs(fftshift(fft(VC1(1:LL)))/LL)));grid on;hold on;xlim([-1e10,1e10]);
% figure(8),plot(f,log10(abs(fftshift(fft(VC2(1:LL)))/LL)));grid on;hold on;xlim([-1e10,1e10]);
% figure(9),plot(VC1(1:LL));grid on;hold on;
% clear