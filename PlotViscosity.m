% viscosity, different means

clear
clf

% EDIT THESE THINGS
directory = 'D:\WSLFiles\Output\Venus\210507_testStagnentLid\+op';
file_name = 'Venus';
file_stem = [directory file_name];
plate_analyse = strcat(directory, '\',file_name, '_plates_analyse.dat');

starting_frame = 0;
num = 5; % plot every xxx frames
total_frame = size(importdata(plate_analyse).data, 1)-1
frames = starting_frame:num:total_frame;

viscosity.geometric = zeros(1,length(starting_frame:num:total_frame));
viscosity.arithmetic =  zeros(1,length(starting_frame:num:total_frame));

% Dimensional scales
Dscale=0.001;     % depth scale: m to km
tscale=1/(3600*24*365.24*1e6);    % age scale: s to Myr
Depth = 2942;    % in km

% Load spatial fields: 4D arrays (x,y,z,b)
for i = 1:length(frames)
    frame = frames(i);
    disp(strcat('Frame ', num2str(frame)))
    framestring = numstring5(frame);

    [theta phi z visco] = ReadStag3Dpjt(directory, file_name, frame, 'viscosity'); 
    
    viscosity.geometric(i) = geomean(reshape(squeeze(visco), 1,[]));
    viscosity.arithmetic(i) = mean(reshape(squeeze(visco), 1, []));
end

semilogy(frames,viscosity.geometric,frames,viscosity.arithmetic)
legend('geometric mean', 'arithmetic mean','location','southeast')
xlabel('Frame')
ylabel('Viscosity')