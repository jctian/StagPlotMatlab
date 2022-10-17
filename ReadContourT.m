% Read 2D average temperature/surface age/crust thickness from StagYY outputs
% only 2D case for venus 

clear

start_dir = 'D:\WSLFiles\StagPlotMatlab';
addpath('D:\WSLFiles\StagPlotMatlab')
directory = 'D:\Output\Venus_2022\220708_noharzmelt_erupt100ol\+op';
file_name = 'Venus';
nx = 512;
nz = 96;

plot_fig = true;

try
    plate_analyse = strcat(directory, '\',file_name, '_plates_analyse.dat');
    time = strcat(directory, '\',file_name, '_time.dat');
    plate_analyse_data = importdata(plate_analyse);
    time_data = importdata(time);
    number_frame = size(plate_analyse_data.data, 1)-3; % in output file should be step 0 to number_frame-1
    time_full = time_data.data(:,2); % first colume frame, second time
    step_full = time_data.data(:,1);
    
    refstate = strcat(directory, '\',file_name, '_refstat.dat');
    combinedadiabat = zeros(192,8); % z    T     rho     expan      Cp     tcond     P   grav
    
    % % number of lines for refstate
    % fid = fopen(refstate, 'rt');
    % chunksize = 1e6; % read chuncks of 1MB at a time
    % numRows = 1;
    % while ~feof(fid)
    %    ch = fread(fid, chunksize, '*uchar');
    %    if isempty(ch)
    %        break
    %    end
    %    numRows = numRows + sum(ch == sprintf('\n'));
    % end
    % fclose(fid);
    
    % skip until last 96+1 lines (empty line at last)
    % in = 1;
    % fid2 = fopen(refstate, 'rt');
    % for i = 1:numRows-1
    %     if i < numRows-193
    %         fgetl(fid2);
    %     else
    %         str = fgetl(fid2);
    %         combinedadiabat(in,:) = str2num(str);
    %         in = in+1;
    %     end
    % end
    
    
    frames = plate_analyse_data.data(:,1); % steps for frames
    time_frame = zeros(number_frame,1) ; % time when each frame is recorded
    
    for i = 1:number_frame
        for j = 1:length(step_full)
            if step_full(j) == frames(i)
                time_frame(i) = time_full(j)/(3600*24*365.24*1e9); % to Ga
            end
        end
    end 
    
    % if(step_full(length(step_full))+1 == length(step_full)) % every frame is recorded in the file
    %     for i = 1:number_frame
    %         time_frame(i) = time_full(frames(i)+1)/(3600*24*365.24*1e9); % to Ga
    %         
    %     end
    % else
    %     error("check frames")
    % end
    
    age_ts = zeros(1,number_frame);    % ts = time series
    crthick_ts = zeros(1,number_frame);
    T_ts = zeros(nz,number_frame);  % average temperature, in radical (so 1*nz*totalTimeStep)
    BS_ts = zeros(nz,number_frame); % average basalt
    sedot_ts = zeros(nx,number_frame); % average surface strain
    meltcr_ts = zeros(nx, number_frame); % average melt above the grid where the base of the crust is 
    age_theta_ts = zeros(nx,number_frame); %surface age 
    
    UMLM_boundary = 730e3;  % which phase change? for ol or px-gt
    mantle_D = 2942e3;     % mantle depth/thickness
    
    % Dimensional scales
    Dscale=0.001;     % depth scale: m to km
    tscale=1/(3600*24*365.24*1e6);    % age scale: s to Myr
    vscale=100*3600*24*365.24;    % velocity scale: m/s to cm/year
    
    step = 1;
    
    if strcmp(pwd, directory) == 0
        start_dir = pwd;
    else
        start_dir = 'D:\Polybox\ETH\2020Fall\Research\Plotting\StagVTKMatlab';
    end
    cd(directory) % time-comsuming, should be done in main file (or modify readstag file and add folder directory to all file paths)
    
    for frame = 0:step:(number_frame-1)
        disp(['Frame ', num2str(frame)])
        % reading surface age
        [theta phi z age] = ReadStag3Dpjt(directory, file_name, frame, 'age'); % all z levels
        nz    = size(age,3); 
        sage  = squeeze(age(:,:,nz,:))*tscale;  % age of outermost level ('surface')
        nhp   = prod(size(sage));
        age1d = squeeze(reshape(sage ,nhp,1));  % 1-D scaled array for making histograms
        age_ts(frame+1)=mean(age1d);
        age_theta_ts(:,frame+1) = sage;
        
        %reading crust thickness
        [theta phi z crdat] = ReadStag3Dpjt(directory, file_name, frame, 'crustal thickness'   ); % has nz=1
        crthick= squeeze(crdat)*Dscale;  % get rid of nz=1 dimension
        nhp    = prod(size(crthick));
        cr1d   = squeeze(reshape(crthick ,nhp,1));  % convert to 1-D array
        crthick_ts(frame+1) = mean(cr1d);
        
        %reading temperature
        [theta phi z temperature] = ReadStag3Dpjt(directory, file_name, frame, 'temperature'); % all z levels
        stemp = squeeze(temperature);
        stemp_mean = mean(stemp,1); % mean of each layer 
        T_ts(:,frame+1) = stemp_mean'; 
        
        % reading basalt
        [theta phi z basalt] = ReadStag3Dpjt(directory, file_name, frame, 'basalt'); % all z levels
        sbasalt = squeeze(basalt);    
        sbasalt_mean = mean(sbasalt,1); % mean of each layer 
        BS_ts(:,frame+1) = sbasalt_mean'; 
        
        % read strain rate
        [theta phi z edot] = ReadStag3Dpjt(directory, file_name, frame, 'strain rate'         ); % all z levels
        nz     = size(edot,3); 
        sedot  = squeeze(edot(:,:,nz,:));  % strain rate of outermost level ('surface')
        sedot_ts(:,frame+1)=sedot';
        
        % read melting in crust
        [theta phi z meltfrac] = ReadStag3Dpjt(directory, file_name, frame, 'melt fraction'                ); % all z levels
        nz    = size(meltfrac,3); 
        smelt = zeros(size(meltfrac(:,:,nz,:)));
        for j = 1:size(meltfrac,2)
            for i = 1:size(meltfrac,3)
               if (mantle_D-z(i) <= crthick(j)*1000)
                   break
               end
            end
            if i < size(meltfrac,3)
                smelt(j) = max(squeeze(meltfrac(:,j,(i+1):nz,:)))';
            end
        end
        meltcr_ts(:,frame+1) = smelt';
    end
    
    if plot_fig
        figure(1)
        yyaxis left
        plot(time_frame(1:step:number_frame), age_ts(1:step:number_frame))
        
        xlabel('Time/Ga')
        ylabel('Age/Ma')
        xlim([0, 4.5])
        yyaxis right
        plot(time_frame(1:step:number_frame), crthick_ts(1:step:number_frame))
        title('Evolution of Average Surface Age and Crustal Thickness')
        ylabel('Thickness/km')
        
        legend('Surface Age', 'Crustal Thickness')
        xlim([0, 4.5])
        
        mean(cr1d)
        std(cr1d)
        
        figure(2)
        subplot(2,1,1)
        depth = mantle_D/1000 - z/1000;
        [X,Y] = meshgrid(time_frame(1:step:number_frame), depth);
        [C,h] = contourf(X,Y,T_ts,[700:10:4000]);
        
        set(h,'LineColor','none')
        axis ij
        load("vik.mat")
        colormap(vik)
        colorbar
        % colorbar('Ticks',[1000:300:4000])
        % colormap(brewermap(128,'*spectral'))
        xlabel('Time/Ga')
        ylabel('Depth/km')
        title('Average Temperature')
        
        subplot(2,1,2)
        [C,h] = contourf(X,Y,BS_ts,[0:0.01:1]);
        set(h,'LineColor','none')
        axis ij
        % colorbar('Ticks',[0:0.2:1])
        % colormap(brewermap(128,'*spectral'))
        colormap(vik)
        colorbar
        xlabel('Time/Ga')
        ylabel('Depth/km')
        title('Average Basalt Fraction')
        cd(start_dir)
        
        disp('Finish figure 2')
        
        figure(3)
        subplot(3,1,1)
        nphi   = length(phi);
        dph=360/nphi  ; 
        phip  =0.5*dph:dph:360-0.5*dph;
        sedot_log_ts = log10(sedot_ts);
        sedot_max = max(sedot_log_ts, [], 'all');
        sedot_min = min(sedot_log_ts, [], 'all');
        [X,Y] = meshgrid(time_frame(1:step:number_frame), phip);
        [C,h] = contourf(X,Y,sedot_log_ts,[sedot_min:0.1:sedot_max]);
        set(h,'LineColor','none')
        colorbar
        colormap(brewermap(128,'*spectral'))
        % colormap(vik)
        % colorbar
        xlabel('Time/Ga')
        ylabel('Azimuth/Degree')
        title('Surface Strain Rate (1/s) in Log Scale')
        
        disp('Finish figure 3.1')
        
        subplot(3,1,2)
        meltcr_ts(meltcr_ts < 0.1) = 0;
        [C,h] = contourf(X,Y,meltcr_ts, [0:0.1:1]);
        colormap(brewermap(128,'*spectral'))
        colorbar('Ticks',[0:0.2:1])
        % colormap(vik)
        % colorbar
        xlabel('Time/Ga')
        ylabel('Azimuth/Degree')
        title('Melt Fraction (max melt fraction > 0.1, at this location within the crust)')
        disp('Finish figure 3.2')
        
        subplot(3,1,3)
        
        agetheta_max = max(age_theta_ts,[],'all');
        [C,h] = contourf(X,Y,age_theta_ts,[0:10:agetheta_max]);
        set(h,'LineColor','none')
        colorbar
        colormap(brewermap(128,'*spectral'))
        % colormap(vik)
        % colorbar
        xlabel('Time/Ga')
        ylabel('Azimuth/Degree')
        title('Surface Age/Ma')
    end

    rmpath(start_dir)

catch ME
    cd(start_dir)
    rmpath(start_dir)
    rethrow(ME)
end
