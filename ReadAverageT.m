% Read average temperature/surface age/crust thickness from StagYY outputs
% only 2D case

clear

directory = 'E:\StagYYOutput\210620_intrudemidcrust\+op';
file_name = 'Venus';

plate_analyse = strcat(directory, '\',file_name, '_plates_analyse.dat');
plate_analyse_data = importdata(plate_analyse);
total_frame = size(plate_analyse_data.data, 1)-1
time = plate_analyse_data.data(:,2)/(365.25*24*3600*1e9);

age_ts = zeros(1,total_frame);    % ts = time series
crthick_ts = zeros(1,total_frame);
T_cr_ts = zeros(1,total_frame);   % crust temperature
T_um_ts = zeros(1,total_frame);   % upper mantle T, including crust
T_lm_ts = zeros(1,total_frame);   % lower mantle T
T_ts = zeros(1,total_frame);

UMLM_boundary = 730e3;  % which phase change? for ol or px-gt
mantle_D = 2942e3;     % mantle depth/thickness

% Dimensional scales
Dscale=0.001;     % depth scale: m to km
tscale=1/(3600*24*365.24*1e6);    % age scale: s to Myr
vscale=100*3600*24*365.24;    % velocity scale: m/s to cm/year

step = 1;


start_dir = pwd;
cd(directory) % time-comsuming, should be done in main file

for frame = 1:step:total_frame
    disp(['Frame ', num2str(frame)])
    % reading surface age
    [theta phi z age] = ReadStag3Dpjt(directory, file_name, frame, 'age'); % all z levels
    nz    = size(age,3); 
    sage  = squeeze(age(:,:,nz,:))*tscale;  % age of outermost level ('surface')
    nhp   = prod(size(sage));
    age1d = squeeze(reshape(sage ,nhp,1));  % 1-D scaled array for making histograms
    age_ts(frame)=mean(age1d);
    
    %reading crust thickness
    [theta phi z crdat] = ReadStag3Dpjt(directory, file_name, frame, 'crustal thickness'   ); % has nz=1
    crthick= squeeze(crdat)*Dscale;  % get rid of nz=1 dimension
    nhp    = prod(size(crthick));
    cr1d   = squeeze(reshape(crthick ,nhp,1));  % convert to 1-D array
    crthick_ts(frame) = mean(cr1d);
    
    %reading temperature
    [theta phi z temperature] = ReadStag3Dpjt(directory, file_name, frame, 'temperature'); % all z levels
    stemp = squeeze(temperature);
    
    temp_cr = 0;
    temp_cr_count = 0;
    mantle_grid = 0;
    
    for i = 1:size(stemp,1)  % for every column
        for j = 1:size(stemp,2)
            if(mantle_D-z(j)<=cr1d(j)/Dscale)
                temp_cr = temp_cr + sum(stemp(i,j:96),'all');
                temp_cr_count = temp_cr_count + size(j:96,2);
                break
            end
        end
    end
    
    for j=1:size(stemp,2)
        if (z(j) >= UMLM_boundary)
            mantle_grid =j;
            break
        end 
    end
        
    T_ts(frame) = mean(stemp,'all');
    T_cr_ts(frame) = temp_cr/temp_cr_count;
    T_um_ts(frame) = mean(stemp(:,(mantle_grid):96),'all');
    T_lm_ts(frame) = mean(stemp(:,1:(mantle_grid-1)),'all');
    
    % also melt above crust thickness?
end

subplot(3,1,1)
plot(time(1:step:total_frame), age_ts(1:step:total_frame))
title('Average surface age')
subplot(3,1,2)
plot(time(1:step:total_frame), crthick_ts(1:step:total_frame))
title('Average crust thickness')
subplot(3,1,3)
hold on
plot(time(1:step:total_frame), T_ts(1:step:total_frame))
% plot(time(1:step:total_frame), T_cr_ts(1:step:total_frame))
plot(time(1:step:total_frame), T_um_ts(1:step:total_frame))
plot(time(1:step:total_frame), T_lm_ts(1:step:total_frame))
hold off
legend('T','UM', 'LM')
title('Average temperature')

cd(start_dir)
