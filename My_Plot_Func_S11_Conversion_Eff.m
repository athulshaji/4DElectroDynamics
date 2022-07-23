function My_Plot_Func_S11_Conversion_Eff
clc;close all;clear;
%% Eff vs freq
fid1=fopen('Conversion_Eff_0.95GHz_600_Ohms.txt','r');
fid2=fopen('Conversion_Eff_1.79GHz_600_Ohms.txt','r');
fid3=fopen('Conversion_Eff_2.4GHz_600_Ohms.txt','r');
fgetl(fid1);
fgetl(fid2);
fgetl(fid3);
index=1;
while(~feof(fid1))
    str=fgetl(fid1);
    str=regexp(str,'\t','split');
    pin1(index)=str2num(str{1});
    Eff1(index)=str2num(str{2});
    index=index+1;
end
index=1;
while(~feof(fid2))
    str=fgetl(fid2);
    str=regexp(str,'\t','split');
    pin2(index)=str2num(str{1});
    Eff2(index)=str2num(str{2});
    index=index+1;
end
index=1;
while(~feof(fid3))
    str=fgetl(fid3);
    str=regexp(str,'\t','split');
    pin3(index)=str2num(str{1});
    Eff3(index)=str2num(str{2});
    index=index+1;
end
fclose('all');
figure;
plot(pin1,Eff1,pin2,Eff2,pin3,Eff3,'LineWidth',2);
grid on;
legend('0.95GHz','1.79GHz','2.4GHz');
xlabel('input power');
ylabel('Conversion efficieny');
title('Conversion efficiency vs input power');
axis tight
clear;

%% S11 vs freq (for different load impedances)
fid1=fopen('LSSP_S11_0dBm_200_Ohms.txt','r');
fid2=fopen('LSSP_S11_0dBm_400_Ohms.txt','r');
fid3=fopen('LSSP_S11_0dBm_1000_Ohms.txt','r');
fgetl(fid1);
fgetl(fid2);
fgetl(fid3);
index=1;
while(~feof(fid1))
    str=fgetl(fid1);
    str=regexp(str,'\t','split');
    freq1(index)=str2num(str{1});
    LSSP1(index)=str2num(str{2});
    index=index+1;
end
index=1;
while(~feof(fid2))
    str=fgetl(fid2);
    str=regexp(str,'\t','split');
    freq2(index)=str2num(str{1});
    LSSP2(index)=str2num(str{2});
    index=index+1;
end
index=1;
while(~feof(fid3))
    str=fgetl(fid3);
    str=regexp(str,'\t','split');
    freq3(index)=str2num(str{1});
    LSSP3(index)=str2num(str{2});
    index=index+1;
end
fclose('all');
figure;
plot(freq1,LSSP1,freq2,LSSP2,freq3,LSSP3,'LineWidth',2);
grid on;
legend('200 \Omega','600 \Omega','1000 \Omega');
xlabel('Freq(GHz)');
ylabel('S_{11}(dB)');
title('S_{11}(dB) vs freq');
axis tight
clear;

%% S11 vs freq (for different power input)
fid1=fopen('LSSP_S11_-3dBm_400_Ohms.txt','r');
fid2=fopen('LSSP_S11_0dBm_400_Ohms.txt','r');
fid3=fopen('LSSP_S11_3dBm_400_Ohms.txt','r');
fgetl(fid1);
fgetl(fid2);
fgetl(fid3);
index=1;
while(~feof(fid1))
    str=fgetl(fid1);
    str=regexp(str,'\t','split');
    freq1(index)=str2num(str{1});
    LSSP1(index)=str2num(str{2});
    index=index+1;
end
index=1;
while(~feof(fid2))
    str=fgetl(fid2);
    str=regexp(str,'\t','split');
    freq2(index)=str2num(str{1});
    LSSP2(index)=str2num(str{2});
    index=index+1;
end
index=1;
while(~feof(fid3))
    str=fgetl(fid3);
    str=regexp(str,'\t','split');
    freq3(index)=str2num(str{1});
    LSSP3(index)=str2num(str{2});
    index=index+1;
end
fclose('all');
figure;
plot(freq1,LSSP1,freq2,LSSP2,freq3,LSSP3,'LineWidth',2);
grid on;
legend('-3 dBm','0 dBm','3  dBm');
xlabel('Freq(GHz)');
ylabel('S_{11}(dB)');
title('S_{11}(dB) vs freq');
axis tight
clear;
end