% Read in various binary and ascii files and plot crustal thickness, 
%  surface topography, geoid, self-gravitating surface topography, age, and
%  strain rate. 
% Assumes geometry is a full sphere, either 3D yin-yang or 2D annulus.
% Paul Tackley, February 2021

% EDIT THESE THINGS

function PlotSurfFields_function(directory, file_name, frame, include)
    clf
    
    lw=1;   % line width
    
    file_stem = [directory file_name];

    % Dimensional scales
    Dscale=0.001;     % depth scale: m to km
    tscale=1/(3600*24*365.24*1e6);    % age scale: s to Myr
    Depth = 2942;    % in km


    framestring = numstring5(frame);

    cslmname= strcat(file_stem,'_cslm',framestring,'.dat');
    glmname = strcat(file_stem, '_glm',framestring,'.dat');
    YinYang = false;
    ngraphs=0;

    % Load spatial fields: 4D arrays (x,y,z,b)
    if include.topography
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
    if include.crustthick
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
    if include.age
        [theta phi z age] = ReadStag3Dpjt(directory, file_name, frame, 'age'                 ); % all z levels
        nz    = size(age,3); 
        if (include.crustage && include.crustthick)
            % only for 2D
            sage = zeros(size(age(:,:,nz,:)));
            for j = 1:size(age,2)
                for i = 1:size(age,3)
                   if (Depth-z(i)*Dscale <= crthick(j))
                       break
                   end
                end
                sage(j) = mean(squeeze(age(:,j,i:nz,:))).*tscale';
            end
        else
            sage  = squeeze(    age(:,:,nz,:))*tscale;  % age of outermost level ('surface')
        end
        
        nhp   = prod(size(sage));
        age1d = squeeze(reshape( sage ,nhp,1));  % 1-D scaled array for making histograms
        if(nhp>max(size(sage))) % detect 3D
            YinYang = true;
            agemap= YYtoMap2(sage);
        end
        ngraphs= ngraphs+2;
    end
    if include.strainrate
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
    if include.geoid
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

    if include.surfacemelt
        [theta phi z meltfrac] = ReadStag3Dpjt(directory, file_name, frame, 'melt fraction'                ); % all z levels
        nz    = size(meltfrac,3); 
        if include.crustmelt
            % only for 2D
            smelt = zeros(size(meltfrac(:,:,nz,:)));
            smelt_1 = smelt;
            for j = 1:size(meltfrac,2)
                for i = 1:size(meltfrac,3)
                   if (Depth-z(i)*Dscale <= crthick(j))
                       break
                   end
                end
                smelt(j) = mean(squeeze(meltfrac(:,j,i:nz,:)))';  % i+1?
                smelt_1(j) = mean(squeeze(meltfrac(:,j,(i+1):nz,:)))';
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
    
    
%     
%     [theta phi z VX VY VZ P time rcmb] = ReadStag3Dpjt(directory, file_name, frame, 'velocity');
%     VX = squeeze(VX);
%     VY = squeeze(VY);
%     VX_surf = VX(:,90:96);
%     VY_surf = VY(:,90:96);

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

    if include.topography
        graph = graph+1;
        subplot(nrows,2,graph); 
        if YinYang
            contourf(phip,thetap,smap); %hold on; plot(theta,sg,'r'); hold off; 
            colorbar; ylabel('theta');
        else
            plot(phip,stopo); ylabel('Elevation (km)');
        end
        title('Surface Topography'); %legend('raw','sg filtered')
        xlabel('phi')
        graph = graph+1;
        subplot(nrows,2,graph);
        hist(topo1d,nbins); title('Surface Topography'); xlabel('elevation (km)')
    end
    if include.age
        graph = graph+1;
        subplot(nrows,2,graph); 
        if YinYang
            contourf(phip,thetap,agemap); 
            colorbar; ylabel('theta');
        else
            plot(phip,sage); ylabel('Age (Myr)');
        end
        if include.crustage
            title('Average crust age');
        else
            title('Surface Age');
        end
        xlabel('phi')
        graph = graph+1;
        subplot(nrows,2,graph);
        hist(age1d,nbins); title('Surface Age'); xlabel('age (Ma)')
    end
    if include.strainrate
        graph = graph+1;
        subplot(nrows,2,graph); 
        if YinYang
            contourf(phip,thetap,edmap); colorbar; ylabel('theta')
        else
            plot(phip,sedot); ylabel('Strain rate (/s)');
        end
        title('Surface Strain Rate'); xlabel('phi')
        graph = graph+1;
        subplot(nrows,2,graph);
        hist(ed1d,nbins); title('Surface Strain Rate'); xlabel('strain rate (/s)')
    end

    if include.surfacemelt
        graph = graph+1;
        subplot(nrows,2,graph); 
        if YinYang
            contourf(phip,thetap,meltmap); 
            colorbar; ylabel('theta');
        else
            plot(phip,smelt,'linewidth',lw); ylabel('Melt fraction');
        end
        if include.crustmelt
            title('Average melt fraction')
            hold on
            plot(phip,smelt_1)
            hold off
            legend('above intrusion depth', 'within the crust')
        else
            title('Surface melt fraction')
        end
        xlabel('phi')
%         graph = graph+1;
%         subplot(nrows,2,graph);
%         hist(melt1d,[0:0.1:1]);
%         if include.crustmelt
%             title('Average melt fraction within crust'); xlabel('Melt fraction')
%         else
%             title('Surface melt fraction'); xlabel('Melt fraction')
%         end
%         xlim([0,1])
    end

    if include.crustthick
        graph = graph+1;
        subplot(nrows,2,graph); 
        if YinYang
            contourf(phip,thetap,crmap); colorbar; ylabel('theta')
        else
            plot(phip,crthick); ylabel('Thickness (km)');
        end
        title('Crustal Thickness'); xlabel('phi')
    end
    if include.geoid
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
    
    sgtitle(strcat('Frame ', num2str(frame))) 
    sgt.FontSize = 12;
    pause(0.5)
    set(gcf,'Position',[100 100 800 800])

end
