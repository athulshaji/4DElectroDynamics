function MyDataVisualizer
clc;
% close all;
clear variables;
% code used to generate the results from the TD3DRad.m code

VIDEO = 0;
if(VIDEO)
    %     fig = figure('Color', 'w', 'units', 'normalized', 'OuterPosition',[0 0 1 1]);
    fig=figure;
    movieName = 'FieldYSurfDiodeSMS7630.avi';
    videoObj = VideoWriter(movieName);
    open(videoObj);
end
tic
% FieldX=load('FieldXSurf.mat');
DATA=load('FieldZDiodeSMS7630_SPICE_LC_1uF_400Ohm_SinusoidalWave_withoutGNDplane.mat');
[pks,locs] = findpeaks(DATA.FieldZ(:,32));
% DATA=load('FieldZDiodeCap_0.01uF_180_90.mat');
delt=load('delt.txt');
T=2^13;
data = DATA.FieldZ(:,32);
Fs=1/delt;Y = fft(data,T);
L=T;
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
[pks,locs] = findpeaks(P1);
amplitudeFreq = [pks,f(locs).'/1e9];
figure
subplot 211
plot(f/1e9,((P1))); grid on; axis tight;
title('FFT Ey SPICE 1uF');xlabel('GHz');ylabel('FFT mag');
subplot 212
plot((0:T-1)*delt, data); grid on; axis tight;
title('Waveform amplitude vs time');
% FieldZ=load('FieldZSurf.mat');
toc
sz=size(DATA.FieldY);
DATA = load('FieldYSurfDiodeSMS7630_0.01uF.mat');
sz = size(DATA.FieldYSurf);
if(VIDEO)
    for i=1:T
        imageData=reshape(DATA.FieldYSurf(i,:,:),sz(2),sz(3));
        pcolor(imageData./max(imageData(:))); shading interp;title('iter #',num2str(i));
        caxis([-1,1]);
        drawnow;
        FRAME = getframe(fig);
        writeVideo(videoObj, FRAME)
    end
    close(videoObj);
end

end