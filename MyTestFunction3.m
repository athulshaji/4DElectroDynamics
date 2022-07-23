function MyTestFunction3
% This function is used along with HarmonicBalanceSolver.m to plot the
% convergence of the magnitude of the coefficients corresponding to DC,
% first harmonic, and higher harmonics
clc;close all;clear all;
format long;
voptim=load('voptim.mat');
L=length(voptim.voptim);
for i=1:L
    absvoptim{i}=abs(voptim.voptim{i});
end
for i=1:L
    magDC(i)=absvoptim{i}(i+1);
    mag1HAR(i)=absvoptim{i}(i+2);
end
meanmagDC=mean(magDC);
stdmagDC=std(magDC);
meanmag1HAR=mean(mag1HAR);
stdmag1HAR=std(mag1HAR);
xlabel_name='# harmonics';
figure('Color','w','units','normalized','outerposition',[0 0 1 1]);
subplot 241,
yyaxis left
plot(magDC,'-ok','LineWidth',2);axis tight;grid on;xlabel(xlabel_name);
ylabel('magnitude');title('DC value');
text(10,max(magDC)-stdmagDC,['mean=',num2str(meanmagDC)]);
text(10,max(magDC)-1.5*stdmagDC,['std deviation=',num2str(stdmagDC)])
yyaxis right
ylabel('deviation from mean');
plot(abs(meanmagDC-magDC),'LineWidth',2);
subplot 242,
yyaxis left
plot(mag1HAR,'-ok','LineWidth',2);axis tight;grid on;xlabel(xlabel_name);
ylabel('magnitude');title('First harmonic');
text(10,max(mag1HAR)-stdmag1HAR,['mean=',num2str(meanmag1HAR)]);
text(10,max(mag1HAR)-1.5*stdmag1HAR,['std deviation=',num2str(stdmag1HAR)]);
yyaxis right
ylabel('deviation from mean');
plot(abs(meanmag1HAR-mag1HAR),'LineWidth',2);
for i=2:L
    mag2HAR(i)=absvoptim{i}(i+3);
end
subplot 243,
yyaxis left
plot(mag2HAR(2:L),'-ok','LineWidth',2);axis tight;grid on;xlabel(xlabel_name);
ylabel('magnitude');title('2nd harmonic');
meanmag2HAR=mean(mag2HAR);
stdmag2HAR=std(mag2HAR);
text(10,max(mag2HAR)-stdmag2HAR,['mean=',num2str(meanmag2HAR)]);
text(10,max(mag2HAR)-1.15*stdmag2HAR,['std deviation=',num2str(stdmag2HAR)]);
yyaxis right
ylabel('deviation from mean');
plot(abs(meanmag2HAR-mag2HAR(2:L)),'LineWidth',2);
for i=3:L
    mag3HAR(i)=absvoptim{i}(i+4);
end
subplot 244,
yyaxis left
plot(mag3HAR(3:L),'-ok','LineWidth',2);axis tight;grid on;xlabel(xlabel_name);
ylabel('magnitude');title('3rd harmonic');
meanmag3HAR=mean(mag3HAR);
stdmag3HAR=std(mag3HAR);
text(10,max(mag3HAR)-0.25*stdmag3HAR,['mean=',num2str(meanmag3HAR)]);
text(10,max(mag3HAR)-0.3*stdmag3HAR,['std deviation=',num2str(stdmag3HAR)]);
yyaxis right
ylabel('deviation from mean');
plot(abs(meanmag3HAR-mag3HAR(3:L)),'LineWidth',2);
for i=4:L
    mag4HAR(i)=absvoptim{i}(i+5);
end
yyaxis left
subplot 245,
yyaxis left
plot(mag4HAR(4:L),'-ok','LineWidth',2);axis tight;grid on;xlabel(xlabel_name);
ylabel('magnitude');title('4th harmonic');
meanmag4HAR=mean(mag4HAR);
stdmag4HAR=std(mag4HAR);
text(10,max(mag4HAR)-0.25*stdmag4HAR,['mean=',num2str(meanmag4HAR)]);
text(10,max(mag4HAR)-0.35*stdmag4HAR,['std deviation=',num2str(stdmag4HAR)]);
yyaxis right
ylabel('deviation from mean');
plot(abs(meanmag4HAR-mag4HAR(4:L)),'LineWidth',2);
for i=5:L
    mag5HAR(i)=absvoptim{i}(i+6);
end
subplot 246,
yyaxis left
plot(mag5HAR(5:L),'-ok','LineWidth',2);axis tight;grid on;xlabel(xlabel_name);
ylabel('magnitude');title('5th harmonic');
meanmag5HAR=mean(mag5HAR);
stdmag5HAR=std(mag5HAR);
text(10,max(mag5HAR)-stdmag5HAR,['mean=',num2str(meanmag5HAR)]);
text(10,max(mag5HAR)-1.5*stdmag5HAR,['std deviation=',num2str(stdmag5HAR)]);
yyaxis right
ylabel('deviation from mean');
plot(abs(meanmag5HAR-mag5HAR(5:L)),'LineWidth',2);
for i=6:L
    mag6HAR(i)=absvoptim{i}(i+7);
end
subplot 247,
yyaxis left
plot(mag6HAR(6:L),'-ok','LineWidth',2);axis tight;grid on;xlabel(xlabel_name);
ylabel('magnitude');title('6th harmonic');
meanmag6HAR=mean(mag6HAR);
stdmag6HAR=std(mag6HAR);
text(10,max(mag6HAR)-0.25*stdmag6HAR,['mean=',num2str(meanmag6HAR)]);
text(10,max(mag6HAR)-0.3*stdmag6HAR,['std deviation=',num2str(stdmag6HAR)]);
yyaxis right
ylabel('deviation from mean');
plot(abs(meanmag6HAR-mag6HAR(6:L)),'LineWidth',2);
for i=7:L
    mag7HAR(i)=absvoptim{i}(i+8);
end
subplot 248,
yyaxis left
plot(mag7HAR(7:L),'-ok','LineWidth',2);axis tight;grid on;xlabel(xlabel_name);
ylabel('magnitude');title('7th harmonic');
meanmag7HAR=mean(mag7HAR);
stdmag7HAR=std(mag7HAR);
text(10,max(mag7HAR)-0.25*stdmag7HAR,['mean=',num2str(meanmag7HAR)]);
text(10,max(mag7HAR)-0.3*stdmag7HAR,['std deviation=',num2str(stdmag7HAR)]);
yyaxis right
ylabel('deviation from mean');
plot(abs(meanmag7HAR-mag7HAR(7:L)),'LineWidth',2);
end