function MyDataVisualizer2
clc;
close all;
clear;
DATAXDiodeCap1=load('FieldXDiodeSMS7630.mat');
DATAYDiodeCap1=load('FieldYDiodeCapSMS7630.mat');
DATAZDiodeCap1=load('FieldZDiodetCapSMS7630.mat');
DATAXDiodeCap2=load('FieldXDiodeCapSMS7630.mat');
DATAYDiodeCap2=load('FieldYDiodeCapSMS7630.mat');
DATAZDiodeCap2=load('FieldZDiodeCapSMS7630.mat');
delt=load('delt.txt');
T=2^13;
positionX = 1:T;
positionY = 32;
FIELD_SHORT_MAG = sqrt(DATAXDiodeCap1.FieldX(positionX,positionY).^2+DATAYShort.FieldY(positionX,positionY).^2+DATAZDiodeCap.FieldZ(positionX,positionY).^2);
FIELD_MAG = sqrt((DATAX.FieldX(positionX,positionY).^2+DATAYDiodeCap2.FieldY(positionX,positionY).^2+DATAZDiodeCap2.FieldZ(positionX,positionY).^2));
FIELD_MAG_NORM = FIELD_MAG./max(FIELD_MAG);
FIELD_SHORT_MAG_NORM = FIELD_SHORT_MAG./max(FIELD_SHORT_MAG);
REL_ERROR = abs((FIELD_SHORT_MAG_NORM - FIELD_MAG_NORM)./FIELD_MAG_NORM);
ABS_ERROR = abs((FIELD_SHORT_MAG_NORM - FIELD_MAG_NORM));
ERROR = ABS_ERROR;
figure('Color','w');
subplot 211
plot((positionX)*delt,ERROR,'LineWidth',3,'Color','k');grid on
title('Abolute error of normalized fields vs time','FontSize',12);
xlabel('Time (sec)','FontSize',12);ylabel('Error','FontSize',12);
subplot 212
plot((positionX)*delt,FIELD_MAG,'LineWidth',3),hold on,plot((positionX)*delt,FIELD_SHORT_MAG,'LineWidth',3);grid on
title('Field mag','FontSize',12);xlabel('Time (sec)','FontSize',12);ylabel('Field mag','FontSize',12);
legend('Field mag for continuous tx line','Field mag for shorted tx line');
end