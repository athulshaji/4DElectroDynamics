function MyFunction_Z
clc;close all;clear;
fid=fopen('inputZ2.s1p');
% fgetl(fid);
% fgetl(fid);
index=1;
while(~feof(fid))
    str=fgetl(fid);
%     if(~strcmp(str,'END'))
        str=regexp(str,' ','split');
        freq(index)=str2num(str{1});
        mreal=str2num(str{2});
        mimag=str2num(str{3});
        realZ(index)=mreal;
        imagZ(index)=mimag;
        index=index+1;
%     end
end
fclose(fid);
data=[freq.' realZ.' imagZ.'];
% save Zin_mod.txt data -ascii
freq=[0 0.3 0.6 0.7 0.8 0.88 0.9 0.925 1.0...
    1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.95 2.0 2.1 2.2 2.3 2.4 2.45 2.5 2.6];
querry=linspace(0,2.6,1001);
realZ=spline(freq,realZ,querry);
imagZ=spline(freq,imagZ,querry);
freq=querry;
figure,subplot 211;plot(freq,realZ);grid on;subplot 212;plot(freq,imagZ);grid on;
% [theta,rho]=cart2pol(realZ,imagZ);
% data=[freq.' realZ.' imagZ.'];
% data=[freq.' rho.' theta.'];
save inputZ.mdf data -ascii
% fid=fopen('inputZ.mdf','w');
% fwrite(fid,data,'double');
% fclose(fid);
end