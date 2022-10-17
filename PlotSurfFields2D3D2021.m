% Read in various binary and ascii files and plot crustal thickness, 
%  surface topography, geoid, self-gravitating surface topography, age, and
%  strain rate. 
% Assumes geometry is a full sphere, either 3D yin-yang or 2D annulus.
% Paul Tackley, February 2021

clear

start_dir = 'D:\Polybox\ETH\2020Fall\Research\Plotting\PlotSurfaceFields_new'
% EDIT THESE THINGS
directory = 'E:\Output\Venus\210714_Oldiff50MPa\210714_Oldiff50MPa_erupt100\+op\';
file_name = 'Venus';

[fields, total_frame] = readNames(directory);
frame = total_frame;


file_stem = [directory file_name];

include_topography=true;
include_age       =true;
include_strainrate=false;
include_crustthick=true;
include_geoid     =false;
include_surfacemelt = false;
include_crustmelt = false; % sum of melt within crust



% Dimensional scales
Dscale=0.001;     % depth scale: m to km
tscale=1/(3600*24*365.24*1e6);    % age scale: s to Myr
Depth = 2942;    % in km


framestring = numstring5(frame);

cslmname= strcat(file_stem,'_cslm',framestring,'.dat');
glmname = strcat(file_stem, '_glm',framestring,'.dat');
YinYang = false;
ngraphs=0;


cd(directory) % time-comsuming, should be done in main file


% Load spatial fields: 4D arrays (x,y,z,b)
if include_topography
    [theta phi z cstopo] = ReadStag3Dpjt(directory, file_name, frame, 'topography'          ); % has nz=2
    stopo  = squeeze(cstopo(:,:, 2,:))*Dscale;  % surface topography is in iz=2
    nhp    = prod(size(stopo));
    topo1d = squeeze(reshape(stopo ,nhp,1));  % convert to 1-D array
    if(nhp>max(size(stopo))) % detect 3D
        YinYang = true;
        smap   = YYtoMap2(stopo);
    end
    ngraphs= ngraphs+2;
end
if include_age
    [theta phi z age] = ReadStag3Dpjt(directory, file_name, frame, 'age'                 ); % all z levels
    nz    = size(age,3); 
    sage  = squeeze(    age(:,:,nz,:))*tscale;  % age of outermost level ('surface')
    nhp   = prod(size(sage));
    age1d = squeeze(reshape( sage ,nhp,1));  % 1-D scaled array for making histograms
    if(nhp>max(size(sage))) % detect 3D
        YinYang = true;
        agemap= YYtoMap2(sage);
    end
    ngraphs= ngraphs+2;
end
if include_strainrate
    [theta phi z edot] = ReadStag3Dpjt(directory, file_name, frame, 'strain rate'         ); % all z levels
    nz     = size(edot,3); 
    sedot  = squeeze(   edot(:,:,nz,:))       ;  % strain rate of outermost level ('surface')
    nhp    = prod(size(sedot));
    ed1d   = squeeze(reshape(sedot ,nhp,1));  % 1-D scaled array for making histograms
    if(nhp>max(size(sedot))) % detect 3D
        YinYang = true;
        edmap  = YYtoMap2(sedot);
    end
    ngraphs= ngraphs+2;
end
if include_crustthick
    [theta phi z crdat] = ReadStag3Dpjt(directory, file_name, frame, 'crustal thickness'   ); % has nz=1
    crthick= squeeze(crdat)*Dscale;  % get rid of nz=1 dimension
    nhp    = prod(size(crthick));
    cr1d   = squeeze(reshape(crthick ,nhp,1));  % convert to 1-D array
    if(nhp>max(size(crthick))) % detect 3D
        YinYang = true;
        crmap  = YYtoMap2(crthick);
    end
    ngraphs= ngraphs+1;
end
if include_geoid
    [theta phi z csgeoid] = ReadStag3Dpjt(directory, file_name, frame, 'geoid'               ); % has nz=2
    [theta phi z csgtopo] = ReadStag3Dpjt(directory, file_name, frame, 'topography self-grav'); % has nz=2
    sgeoid  = squeeze(csgeoid(:,:, 2,:))       ;  % surface geoid too
    sgtopo  = squeeze(csgtopo(:,:, 2,:))*Dscale;  % self-gravitating topography (& filtered)
    nhp     = prod(size(sgeoid));
    gtopo1d = squeeze(reshape(sgtopo,nhp,1));
    if(nhp>max(size(sgeoid))) % detect 3D
        YinYang = true;
        gmap    = YYtoMap2(sgeoid);
        sgmap   = YYtoMap2(sgtopo);
    end
    
    % load spectral fields (topo & geoid)
    %  these files contain lines with l m dcmb_cos dcmb_sin dsurf_cos dsurf_sin
    cslm  = load(cslmname); l = cslm(:,1); lmax = max(l);
    glm   = load(glmname);
    % calculate sum of squared amplitudes for surface topography and geoid
    % need factors of 2 for m>0 ???
    ds_sumsq = zeros(lmax,1);
    dg_sumsq = ds_sumsq;
    dsdg_sum = ds_sumsq;
    for i = 1:length(l)
        if l(i)>0
            ds_sumsq(l(i)) = ds_sumsq(l(i)) + cslm(i,5)^2 + cslm(i,6)^2;
            dg_sumsq(l(i)) = dg_sumsq(l(i)) +  glm(i,5)^2 +  glm(i,6)^2;
            dsdg_sum(l(i)) = dsdg_sum(l(i)) + cslm(i,5)*glm(i,5)+ cslm(i,6)*glm(i,6);
        end
    end
    for i = 1:lmax
        ds_sumsq(i) = ds_sumsq(i) / (2*i+1);
        dg_sumsq(i) = dg_sumsq(i) / (2*i+1);
        dsdg_sum(i) = dsdg_sum(i) / (2*i+1);
    end
    ds_rms = sqrt(ds_sumsq);
    dg_rms = sqrt(dg_sumsq);
    admittance  = dsdg_sum ./ ds_sumsq;
    correlation = dsdg_sum ./ (ds_rms.*dg_rms);
    ngraphs= ngraphs+3;
end

if include_surfacemelt
    [theta phi z meltfrac] = ReadStag3Dpjt(directory, file_name, frame, 'melt fraction'                ); % all z levels
    nz    = size(meltfrac,3); 
    if include_crustmelt
        % only for 2D
        smelt = zeros(size(meltfrac(:,:,nz,:)));
        smelt_1 = smelt;
        for j = 1:size(meltfrac,2)
            for i = 1:size(meltfrac,3)
               if (Depth-z(i)*Dscale <= crthick(j))
                   break
               end
            end
            smelt(j) = sum(squeeze(meltfrac(:,j,i:nz,:)))./i';  % i+1?
            if(i < size(meltfrac,3))
                smelt_1(j) = sum(squeeze(meltfrac(:,j,(i+1):nz,:)))./i';
            end
        end
    else
        smelt  = squeeze(meltfrac(:,:,nz,:));  % melt in surface grid
    end 
    nhp   = prod(size(smelt));
    melt1d = squeeze(reshape(smelt,nhp,1));  % 1-D scaled array for making histograms
    if(nhp>max(size(smelt))) % detect 3D
        YinYang = true;
        meltmap= YYtoMap2(smelt);
    end
    ngraphs= ngraphs+2;
end


% now plot everything

if YinYang
    ntheta = 2*length(theta); 
    nphi   = 4*length(phi)/3; 
else
    ntheta = length(theta);
    nphi   = length(phi);
end
dth=180/ntheta; thetap=0.5*dth:dth:180-0.5*dth;
dph=360/nphi  ; phip  =0.5*dph:dph:360-0.5*dph;

nbins = 20;  % number of bins for histograms
nrows = ceil(ngraphs/2);
graph = 0;

if include_topography
    graph = graph+1;
    subplot(nrows,2,graph); 
    if YinYang
        contourf(phip,thetap,smap); %hold on; plot(theta,sg,'r'); hold off; 
        colorbar; ylabel('theta');
    else
        plot(phip,stopo); ylabel('Elevation/km');
    end
    title('Surface Topography'); %legend('raw','sg filtered')
    xlabel('Azimuth/degrees')
    xlim([0 360])
    graph = graph+1;
    subplot(nrows,2,graph);
    histogram(topo1d,50, 'Normalization', 'probability');
    xlim([-8 8])
    title('Surface Topography'); xlabel('Elevation/km')
end
if include_age
    graph = graph+1;
    subplot(nrows,2,graph); 
    if YinYang
        contourf(phip,thetap,agemap); 
        colorbar; ylabel('theta');
    else
        plot(phip,sage); ylabel('Age/Ma');
    end
    title('Surface Age'); xlabel('Azimuth/degrees')
    xlim([0 360])
    graph = graph+1;
    subplot(nrows,2,graph);
    histogram(age1d,nbins, 'Normalization', 'probability');
    title('Surface Age'); xlabel('Age/Ma')
end
if include_strainrate
    graph = graph+1;
    subplot(nrows,2,graph); 
    if YinYang
        contourf(phip,thetap,edmap); colorbar; ylabel('theta')
    else
        plot(phip,sedot); ylabel('Strain rate (/s)');
    end
    title('Surface Strain Rate');  xlabel('Azimuth/degrees')
    xlim([0 360])
    graph = graph+1;
    subplot(nrows,2,graph);
    hist(ed1d,nbins); title('Surface Strain Rate'); xlabel('Strain Rate (1/s)')
end

if include_surfacemelt
    graph = graph+1;
    subplot(nrows,2,graph); 
    if YinYang
        contourf(phip,thetap,meltmap); 
        colorbar; ylabel('theta');
    else
        plot(phip,smelt); ylabel('Melt fraction');
    end
    if include_crustmelt
        title('Average melt fraction within crust')
        hold on
        plot(phip,smelt_1)
        hold off
        legend('bottom crust', 'bottom crust+1')
    else
        title('Surface melt fraction')
    end
    xlabel('Azimuth/$^\circ$')
    graph = graph+1;
    subplot(nrows,2,graph);
    hist(melt1d,[0:0.1:1]); 
    title('Surface melt fraction'); xlabel('Melt fraction')
    xlim([0,1])
end

if include_crustthick
    graph = graph+1;
    subplot(nrows,2,graph); 
    if YinYang
        contourf(phip,thetap,crmap); colorbar; ylabel('theta')
    else
        plot(phip,crthick); ylabel('Thickness/km');
    end
    title('Crustal Thickness');  xlabel('Azimuth/degrees')
    xlim([0 360])
    graph = graph+1;
    subplot(nrows,2,graph);
    histogram(crthick,nbins, 'Normalization', 'probability'); 
    title('Crustal Thickness'); xlabel('Crustal Thickness/km')
end
if include_geoid
    graph = graph+1;
    subplot(nrows,2,graph); 
    if YinYang
        contourf(phip,thetap,gmap); colorbar; ylabel('theta') 
    else
        plot(phip,sgeoid); ylabel('Geoid (m)');
    end
    title('Geoid'); xlabel('phi')
    graph = graph+1;
    subplot(nrows,2,graph);
      semilogy(ds_rms*Dscale,'r'); hold on; 
      semilogy(dg_rms*Dscale,'g'); hold off; legend('topo','geoid')
      title('Spectra: Topography & Geoid'); ylabel('amplitude (km)'); xlabel('L')
    graph = graph+1;
    subplot(nrows,2,graph)
      plot(correlation,'g'); hold on; plot(admittance,'r'); hold off
      title('Correlation and Admittance Ratios'); xlabel('L'); ylabel('value')
      legend('correlation','admittance')
end

cd(start_dir)

disp('age')
mean(age1d)
std(age1d)
disp('crthick')
mean(crthick)
std(crthick)


