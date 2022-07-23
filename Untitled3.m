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
load('02Oct2021\data\nffieldxx_1.txt');
load('02Oct2021\data\nffieldyy_1.txt');
load('02Oct2021\data\nffieldzz_1.txt');
load('02Oct2021\data\fffieldxx_1.txt');
load('02Oct2021\data\fffieldyy_1.txt');
load('02Oct2021\data\fffieldzz_1.txt');
% case 1: with diode
load('02Oct2021\data\VC1_2.txt');
load('02Oct2021\data\nffieldxx_2.txt');
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