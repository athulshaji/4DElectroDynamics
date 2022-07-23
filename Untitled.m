clc;close all;clear;
format short;
fid=fopen('0dBm_V1.txt');
index=1;
voltage1=[];voltage2=[];time1=[];time2=[];
while(~feof(fid))
    str=fgetl(fid);
    m = str2double(split(str));
    time1(index)=m(1);
    voltage1(index)=m(2);
    index=index+1;
end
fclose(fid);
fid=fopen('0dBm_V2.txt');
index=1;
while(~feof(fid))
    str=fgetl(fid);
    m = str2double(split(str));
    time2(index)=m(1);
    voltage2(index)=m(2);
    index=index+1;
end
fclose(fid);
figure(1),plot(time1,voltage1);hold on;
figure(2),plot(time2,voltage2);hold on;
load('case-0dBm-0.308Vpeak_withVoltageSource.mat');
figure(1),plot(0:delt:(n-1)*delt,V1);
figure(2),plot(0:delt:(n-1)*delt,V2);