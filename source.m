function [tau, t0] = source(fmax)
% fid = fopen(file_name);
% str = fgetl(fid);
% param = str2double(regexp(str, ',', 'split'));
% 
% source_x = param(1); % Excitation Source x coordinate
% source_y = param(2); % Excitation Source y coordinate
% source_z = param(3); % Excitation Source z coordinate

tau=sqrt(2.3)/pi/fmax; % Parameter specifies the Gaussian waveform width in both time and frequency domain
%tau=1/(2*fmax);
t0=5*tau; % The amount of shift requires to avoid non-zero amplitudes at zero time.
end