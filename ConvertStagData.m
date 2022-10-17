    % Read/Write StagYY data
%
% Read binary output files of STAG3D and write them into Paraview files
%
% Boris Kaus, 26.2.2008
%
% PJT edits August 2010:
% - spherical: move coordinate transform just before writing, otherwise reset
% - spherical: centre domain around equator, rotate velocity components
% - YinYang: transform velocity components
% - new Annulus mode (replaced previous 'Cylindrical' mode)

% PJT edits July 2015: Added stress and strain rate
% PJT January 2021: Added basalt, HPE, conductivity, dislocation creep & plasticity

% new main branch: ReadStag3D2 need to add a line for T_core


clear
clc
 
%directory       = "D:\WSLFiles\StagYY\+op"; 
directory       =   'E:\Euler\Stag\Earth_isotope\220926_TestAnnaParFile\+op\'
fname           =   'Gulcher';                      % filename

addpath('D:\WSLFiles\Plotting\StagVTKMatlab')

start_dir = 'D:\WSLFiles\Plotting\StagVTKMatlab';
% scale second to Ma: 3.168876461541279e-14


%GridType        =   'YinYang';         % choose one of these geometries
%GridType        =   'Cartesian';
%GridType        =   'Spherical';
GridType        =   'Annulus';



% Specify which fields we want to read/write in addition to T:
WriteViscosity  =   true;                         % Read & write viscosity files?
WriteVelocity   =   false;                         % Read & write velocity files?
WriteTCond      =   false;
WriteComposition=   false; 
WriteBasalt     =   true;
WriteHPE        =   true;
WriteMeltFrac   =   false;
WriteStress     =   false;
WriteStrainRate =   false;
WriteDislCreep  =   false;
WritePlasticity =   false;
WritePhase      =   false;
WriteMeltComp   =   false;
WriteAge = false;
WriteSmNd142 = false;
WriteHfW182 = false;
% phase field?


try

    [fields, min_frame, max_frame] = readNames(directory); % <3 seconds for a folder with 1400 frames and 17 fields
    
    num = 1; % number of frame for name array
    
    number_start = min_frame;
    number_end = max_frame;
    
    cd(directory)
    
    
    for fname_number=number_start:1:number_end
      fname_number
        switch GridType
            case 'YinYang'
    
                % YING YANG GRID
                disp(['Creating a YIN-YANG grid ...'])
    
                [nnb, X_3D, Y_3D, Z_3D,T_3D_1, T_3D_2, time, rcmb]     = ReadStag3D2(directory, fname, fname_number, 'temperature');
                
                if nnb==-999
                    % The file does not exist, and we have finished processing all
                    % data
    
                    % Create a PVD file, which contains all data including
                    % time-information
                    Create_PVD_file(FileNames,fname,directory);
    
                    error(['finished processing all files in directory'])
    
                else
                    if WriteVelocity
                        % Read pressure & velocity information
                        [nnb, X_3D, Y_3D, Z_3D,VX_3D_1, VY_3D_1, VZ_3D_1, P_3D_1, VX_3D_2, VY_3D_2, VZ_3D_2, P_3D_2, time] = ReadStag3D2(directory, fname, fname_number, 'velocity');
                    end
                    
                    if WriteViscosity
                        % Read viscosity information
                        [nnb, X_3D, Y_3D, Z_3D,ETA_3D_1, ETA_3D_2, time, rcmb] = ReadStag3D2(directory, fname, fname_number, 'viscosity');
                    end
                    
                    if WriteComposition
                        % Read viscosity information
                        [nnb, X_3D, Y_3D, Z_3D,C_3D_1, C_3D_2, time, rcmb]     = ReadStag3D2(directory, fname, fname_number, 'composition');
                    end
                    
                    if WriteMeltFrac
                        % Read viscosity information
                        [nnb, X_3D, Y_3D, Z_3D,F_3D_1, F_3D_2, time, rcmb]     = ReadStag3D2(directory, fname, fname_number, 'melt fraction');
                    end
    
                    if WriteStress
                        % Read stress information
                        [nnb, X_3D, Y_3D, Z_3D,S_3D_1, S_3D_2, time, rcmb]     = ReadStag3D2(directory, fname, fname_number, 'stress');
                    end
    
                    if WriteStrainRate
                        % Read strain rate information
                        [nnb, X_3D, Y_3D, Z_3D,E_3D_1, E_3D_2, time, rcmb]     = ReadStag3D2(directory, fname, fname_number, 'strain rate');
                    end
    
                    if WriteBasalt
                        % Read basalt fraction information
                        [nnb,X_3D, Y_3D, Z_3D,BS_3D_1, BS_3D_2, time]          = ReadStag3D2(directory, fname, fname_number, 'basalt');
                    end
    
                    if WriteHPE
                        % Read Heat-Producing Element information
                        [nnb,X_3D, Y_3D, Z_3D,HPE_3D_1,HPE_3D_2, time]         = ReadStag3D2(directory, fname, fname_number, 'HPE');
                    end
    
                    if WriteTCond
                        % Read thermal conductivity information
                        [nnb,X_3D, Y_3D, Z_3D,K_3D_1,K_3D_2, time]             = ReadStag3D2(directory, fname, fname_number, 'conductivity');
                    end
    
                    if WriteDislCreep
                        % Read disl. creep information
                        [nnb,X_3D, Y_3D, Z_3D,RD_3D_1,RD_3D_2, time]           = ReadStag3D2(directory, fname, fname_number, 'dislocation');
                    end
                    if WritePlasticity
                        % Read plasticity information
                        [nnb,X_3D, Y_3D, Z_3D,RP_3D_1,RP_3D_2, time]           = ReadStag3D2(directory, fname, fname_number, 'plasticity');
                    end
                    
                    xs = X_3D(:,:,end)+pi/4;
                    ys = Y_3D(:,:,end)-pi/4;
                    zs = Z_3D(:,:,end);
    
                    % Transform coordinates for Yin & Yang grids
                    R  = Z_3D+rcmb;     lat = pi/4-X_3D;     lon = Y_3D-3*pi/4;
                    
                    X_3D_1 = R.*cos(lat).*cos(lon);   % Yin grid
                    Y_3D_1 = R.*cos(lat).*sin(lon);
                    Z_3D_1 = R.*sin(lat);
    
                    X_3D_2 = -X_3D_1;       %  Yang grid
                    Y_3D_2 =  Z_3D_1;
                    Z_3D_2 =  Y_3D_1;
                    
                    % Transform velocities, if needed
                    if WriteVelocity
                        Vtheta = VX_3D_1; Vphi = VY_3D_1; Vr = VZ_3D_1;                                 % on Yin grid
                        VX_3D_1 =  Vtheta.*sin(lat).*cos(lon) - Vphi.*sin(lon) + Vr.*cos(lat).*cos(lon);
                        VY_3D_1 =  Vtheta.*sin(lat).*sin(lon) + Vphi.*cos(lon) + Vr.*cos(lat).*sin(lon);
                        VZ_3D_1 = -Vtheta.*cos(lat)                            + Vr.*sin(lat);
                        Vtheta = VX_3D_2; Vphi = VY_3D_2; Vr = VZ_3D_2;                                 % on Yang grid
                        VX_3D_2 = -( Vtheta.*sin(lat).*cos(lon) - Vphi.*sin(lon) + Vr.*cos(lat).*cos(lon) );
                        VZ_3D_2 =  Vtheta.*sin(lat).*sin(lon) + Vphi.*cos(lon) + Vr.*cos(lat).*sin(lon);
                        VY_3D_2 = -Vtheta.*cos(lat)                            + Vr.*sin(lat);
                    end
    
                    % Write a ying-yang grid
                    WriteStag3D_VTK_YinYang2;    % Write VTK XML file
    
                    % Store the filename and time in a structure, which is later
                    % used to create a pvd file, that contains all timesteps and
                    % times of the files
                    FileNames{num} = {time,fname_number,fname_vtk};
    
                end
    
    
            case 'Cartesian'
                % CARTESIAN GRID
                
                % Read temperature information
                [nnb, X_3D, Y_3D, Z_3D,T_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'temperature');
    
                if nnb==-999
                    % The file does not exist, and we have finished processing all
                    % data
    
                    % Create a PVD file, which contains all data including
                    % time-information
                    Create_PVD_file(FileNames,fname,directory);
    
                    error(['finished processing all files in directory'])
    
                else
    
                    if WriteVelocity
                        % Read pressure & velocity information
                        [nnb,X_3D, Y_3D, Z_3D, VX_3D, VY_3D, VZ_3D, P_3D, time] = ReadStag3D2(directory, fname, fname_number, 'velocity');
                    end
    
                    if WriteViscosity
                        % Read viscosity information
                        [nnb,X_3D, Y_3D, Z_3D,ETA_3D, time]   = ReadStag3D2(directory, fname, fname_number, 'viscosity');
                    end
                    
                    if WriteComposition
                        % Read viscosity information
                        [nnb,X_3D, Y_3D, Z_3D,C_3D, time]   = ReadStag3D2(directory, fname, fname_number, 'composition');
                    end
    
                    if WriteMeltFrac
                        % Read viscosity information
                        [nnb,X_3D, Y_3D, Z_3D,F_3D, time]   = ReadStag3D2(directory, fname, fname_number, 'melt fraction');
                    end
                    
                    if WriteStress
                        % Read stress information
                        [nnb,X_3D, Y_3D, Z_3D,S_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'stress');
                    end
    
                    if WriteStrainRate
                        % Read strain rate information
                        [nnb,X_3D, Y_3D, Z_3D,E_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'strain rate');
                    end
                    
                    if WriteBasalt
                        % Read basalt fraction information
                        [nnb,X_3D, Y_3D, Z_3D,BS_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'basalt');
                    end
    
                    if WriteHPE
                        % Read Heat-Producing Element information
                        [nnb,X_3D, Y_3D, Z_3D,HPE_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'HPE');
                    end
    
                    if WriteTCond
                        % Read thermal conductivity information
                        [nnb,X_3D, Y_3D, Z_3D,K_3D, time]      = ReadStag3D2(directory, fname, fname_number, 'conductivity');
                    end
    
                    if WriteDislCreep
                        % Read disl. creep information
                        [nnb,X_3D, Y_3D, Z_3D,RD_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'dislocation');
                    end
                    if WritePlasticity
                        % Read plasticity information
                        [nnb,X_3D, Y_3D, Z_3D,RP_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'plasticity');
                    end    
                    
                    %WriteStag3D_LegacyVTK_Cartesian;   %-> Writes VTK legacy files
                    WriteStag3D_VTK_Cartesian;          %-> Writes VTK-XML files
    
                    % Store the filename and time in a structure, which is later
                    % used to create a pvd file, that contains all timesteps and
                    % times of the files
                    FileNames{num} = {time,fname_number,fname_vtk};
    
                end
                
            case 'Spherical'
                % SPHERICAL GRID
                
                % Read temperature information
                [nnb, X_3D, Y_3D, Z_3D,T_3D, time, rcmb]     = ReadStag3D2(directory, fname, fname_number, 'temperature');
    
                if nnb==-999
                    % The file does not exist, and we have finished processing all
                    % data
    
                    % Create a PVD file, which contains all data including
                    % time-information
                    Create_PVD_file(FileNames,fname,directory);
    
                    error(['finished processing all files in directory'])
    
                else
    
                   if WriteVelocity
                        % Read pressure & velocity information
                        [nnb,X_3D, Y_3D, Z_3D, VX_3D, VY_3D, VZ_3D, P_3D, time] = ReadStag3D2(directory, fname, fname_number, 'velocity');
                    end
    
                    if WriteViscosity
                        % Read viscosity information
                        [nnb,X_3D, Y_3D, Z_3D,ETA_3D, time]   = ReadStag3D2(directory, fname, fname_number, 'viscosity');
                    end
    
                    if WriteComposition
                        % Read viscosity information
                        [nnb,X_3D, Y_3D, Z_3D,C_3D, time]   = ReadStag3D2(directory, fname, fname_number, 'composition');
                    end
    
                    if WriteMeltFrac
                        % Read viscosity information
                        [nnb,X_3D, Y_3D, Z_3D,F_3D, time]   = ReadStag3D2(directory, fname, fname_number, 'melt fraction');
                    end
                    
                    if WriteStress
                        % Read stress information
                        [nnb,X_3D, Y_3D, Z_3D,S_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'stress');
                    end
    
                    if WriteStrainRate
                        % Read strain rate information
                        [nnb,X_3D, Y_3D, Z_3D,E_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'strain rate');
                    end
    
                    if WriteBasalt
                        % Read basalt fraction information
                        [nnb,X_3D, Y_3D, Z_3D,BS_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'basalt');
                    end
    
                    if WriteHPE
                        % Read Heat-Producing Element information
                        [nnb,X_3D, Y_3D, Z_3D,HPE_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'HPE');
                    end
    
                    if WriteTCond
                        % Read thermal conductivity information
                        [nnb,X_3D, Y_3D, Z_3D,K_3D, time]      = ReadStag3D2(directory, fname, fname_number, 'conductivity');
                    end
    
                    if WriteDislCreep
                        % Read disl. creep information
                        [nnb,X_3D, Y_3D, Z_3D,RD_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'dislocation');
                    end
                    if WritePlasticity
                        % Read plasticity information
                        [nnb,X_3D, Y_3D, Z_3D,RP_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'plasticity');
                    end
    
                    % Transform coordinates for Yin grid
                    latmean = (max(max(max(X_3D))) + min(min(min(X_3D)))) / 2;       % spherical patch is centred at equator
                    lonmean = (max(max(max(Y_3D))) + min(min(min(Y_3D)))) / 2;       
                    R    = Z_3D + rcmb;  lat = latmean - X_3D;  lon = Y_3D - lonmean;
                    X_3D = R.*cos(lat).*cos(lon);
                    Y_3D = R.*cos(lat).*sin(lon);
                    Z_3D = R.*sin(lat);
                    
                    % Transform velocities if they are needed
                    if WriteVelocity
                        Vtheta = VX_3D; Vphi = VY_3D; Vr = VZ_3D;
                        VX_3D =  Vtheta.*sin(lat).*cos(lon) - Vphi.*sin(lon) + Vr.*cos(lat).*cos(lon);
                        VY_3D =  Vtheta.*sin(lat).*sin(lon) + Vphi.*cos(lon) + Vr.*cos(lat).*sin(lon);
                        VZ_3D = -Vtheta.*cos(lat)                            + Vr.*sin(lat);
                    end
    
                    %WriteStag3D_LegacyVTK_Cartesian;   %-> Writes VTK legacy files
                    WriteStag3D_VTK_Cartesian;          %-> Writes VTK-XML files
    
                    % Store the filename and time in a structure, which is later
                    % used to create a pvd file, that contains all timesteps and
                    % times of the files
                    FileNames{num} = {time,fname_number,fname_vtk};
    
                end
                
            case 'Annulus'
                % SPHERICAL ANNULUS GRID
                % same as spherical but add extra points to avoid a gap in the
                % visualisation
    
                % Read temperature information
                [nnb, X_3D, Y_3D, Z_3D,T_3D, time, rcmb]     = ReadStag3D2(directory, fname, fname_number, 'temperature');
       
                
                % add one row for nicer visualization
                T_3D(1,end+1,:) = T_3D(1,1,:);
    
                if nnb==-999
                    % The file does not exist, and we have finished processing all
                    % data
    
                    % Create a PVD file, which contains all data including
                    % time-information
                    Create_PVD_file(FileNames,fname,directory);
                    cd(start_dir)
    
                    error(['finished processing all files in directory'])
    
                else
                    
    
                   if WriteVelocity
                        % Read pressure & velocity information
                        [nnb,X_3D, Y_3D, Z_3D, VX_3D, VY_3D, VZ_3D, P_3D, time] = ReadStag3D2(directory, fname, fname_number, 'velocity');
                        VX_3D(1,end+1,:) = VX_3D(1,1,:);
                        VY_3D(1,end+1,:) = VY_3D(1,1,:);
                        VZ_3D(1,end+1,:) = VZ_3D(1,1,:);
                         P_3D(1,end+1,:) =  P_3D(1,1,:);
                   end
    
                    if WriteViscosity
                        % Read viscosity information
                        [nnb,X_3D, Y_3D, Z_3D,ETA_3D, time]   = ReadStag3D2(directory, fname, fname_number, 'viscosity');
                        ETA_3D(1,end+1,:) = ETA_3D(1,1,:);
                    end
    
                    if WriteComposition
                        % Read viscosity information
                        [nnb,X_3D, Y_3D, Z_3D,C_3D, time]   = ReadStag3D2(directory, fname, fname_number, 'composition');
                        C_3D(1,end+1,:) = C_3D(1,1,:);
                    end
    
                    if WriteMeltFrac
                        % Read viscosity information
                        [nnb,X_3D, Y_3D, Z_3D,F_3D, time]   = ReadStag3D2(directory, fname, fname_number, 'melt fraction');
                        F_3D(1,end+1,:) = F_3D(1,1,:);
                    end
                    
                    if WriteStress
                        % Read stress information
                        [nnb,X_3D, Y_3D, Z_3D,S_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'stress');
                        S_3D(1,end+1,:) = S_3D(1,1,:);
                    end
    
                    if WriteStrainRate
                        % Read strain rate information
                        [nnb,X_3D, Y_3D, Z_3D,E_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'strain rate');
                        E_3D(1,end+1,:) = E_3D(1,1,:);
                    end
    
                    if WriteBasalt
                        % Read basalt fraction information
                        [nnb,X_3D, Y_3D, Z_3D,BS_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'basalt');
                        BS_3D(1,end+1,:) = BS_3D(1,1,:);
                    end
    
                    if WriteHPE
                        % Read Heat-Producing Element information
                        [nnb,X_3D, Y_3D, Z_3D,HPE_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'HPE');
                        HPE_3D(1,end+1,:) = HPE_3D(1,1,:);
                    end
    
                    if WriteTCond
                        % Read thermal conductivity information
                        [nnb,X_3D, Y_3D, Z_3D,K_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'conductivity');
                        K_3D(1,end+1,:) = K_3D(1,1,:);
                    end
    
                    if WriteDislCreep
                        % Read disl. creep information
                        [nnb,X_3D, Y_3D, Z_3D,RD_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'dislocation');
                        RD_3D(1,end+1,:) = RD_3D(1,1,:);
                    end
                    if WritePlasticity
                        % Read plasticity information
                        [nnb,X_3D, Y_3D, Z_3D,RP_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'plasticity');
                        RP_3D(1,end+1,:) = RP_3D(1,1,:);
                    end
                    
                    if WritePhase
                        [nnb,X_3D, Y_3D, Z_3D,PH_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'phasefield');
                        PH_3D(1,end+1,:) = PH_3D(1,1,:);
                    end
    
                    if WriteMeltComp
                        [nnb,X_3D, Y_3D, Z_3D,FC_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'melt composition');
                        FC_3D(1,end+1,:) = FC_3D(1,1,:);
                    end
    
                    if WriteAge
                        [nnb,X_3D, Y_3D, Z_3D,AGE_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'age');
                        AGE_3D(1,end+1,:) = AGE_3D(1,1,:);
                    end

                    if WriteSmNd142
                        [nnb,X_3D, Y_3D, Z_3D,SM142_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'Sm142');
                        SM142_3D(1,end+1,:) = SM142_3D(1,1,:);
                        [nnb,X_3D, Y_3D, Z_3D,ND142_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'Nd142');
                        ND142_3D(1,end+1,:) = ND142_3D(1,1,:);
                        [nnb,X_3D, Y_3D, Z_3D,ND144_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'Nd144');
                        ND144_3D(1,end+1,:) = ND144_3D(1,1,:);
                        WriteSm142 = true;
                        WriteNd142 = true;
                        WriteNd144 = true;
                    end
                    
                    if WriteHfW182
                        [nnb,X_3D, Y_3D, Z_3D,W182_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'W182');
                        W182_3D(1,end+1,:) = W182_3D(1,1,:);
                        [nnb,X_3D, Y_3D, Z_3D,HF182_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'Hf182');
                        HF182_3D(1,end+1,:) = HF182_3D(1,1,:);
                        [nnb,X_3D, Y_3D, Z_3D,W184_3D, time]     = ReadStag3D2(directory, fname, fname_number, 'W184');
                        W184_3D(1,end+1,:) = W184_3D(1,1,:);
                        WriteHf182 = true;
                        WriteW182 = true;
                        WriteW184 = true;
                    end
                    
    
    
                    Y_3D(1,end+1,:) = Y_3D(1,1,:) + 2*pi;
                    X_3D(1,end+1,:) = X_3D(1,1,:);
                    Z_3D(1,end+1,:) = Z_3D(1,1,:);
                    
                    % Transform coordinates for Yin grid
                    R    = Z_3D + rcmb;  lat = zeros(size(X_3D));  lon = Y_3D - pi;
                    X_3D = R.*cos(lat).*cos(lon);
                    Y_3D = R.*cos(lat).*sin(lon);
                    Z_3D = R.*sin(lat);
                    
                    % Transform velocities if they are needed
                    if WriteVelocity
                        Vtheta = VX_3D; Vphi = VY_3D; Vr = VZ_3D;
                        VX_3D =  Vtheta.*sin(lat).*cos(lon) - Vphi.*sin(lon) + Vr.*cos(lat).*cos(lon);
                        VY_3D =  Vtheta.*sin(lat).*sin(lon) + Vphi.*cos(lon) + Vr.*cos(lat).*sin(lon);
                        VZ_3D = -Vtheta.*cos(lat)                            + Vr.*sin(lat);
                    end
    
                    %WriteStag3D_LegacyVTK_Cartesian;   %-> Writes VTK legacy files
                    WriteStag3D_VTK_Cartesian;          %-> Writes VTK-XML files
                    
                    if(time == -999)
                        error(['please check input fields'])
                    end
                        
                    % Store the filename and time in a structure, which is later
                    % used to create a pvd file, that contains all timesteps and
                    % times of the files
                    FileNames{num} = {time,fname_number,fname_vtk};
                end
    
            otherwise
    
                error('Unknown gridtype')
        end
    
        num = num+1;
    
    end
    
    Create_PVD_file(FileNames,fname,directory);
    
    cd(start_dir)
    disp('Finished')

catch ME
    cd(start_dir)
    rethrow(ME)
end











