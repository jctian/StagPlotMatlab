% Write STAG_3D format from matlab in XML-VTK format, using ASCII or BINARY file
% format
%
%
%
% Boris Kaus, 21 Feb. 2008
% updated Paul Tackley, July 2015 & January 2021


ASCII           = logical(0);                       % ASCII or BINARY?

% Change to correct directory
%start_dir       = pwd;
cd(directory);


% Generate input data in correct format
Points          = single([X_3D(:) Y_3D(:) Z_3D(:)]);
T_3D            = single(T_3D);


% Transform the arrays in vector-format, and change them to the appropriate
% precision. This example file assumes that data comes in single precision
%
% Double precision requires changing in ths code below (Float32 into Float64)
T               = single(T_3D(:));                                                  % Temperature
Points          = [single(X_3D(:))  single(Y_3D(:))     single(Z_3D(:))     ];      % Coordinates

% Define which arrays to write
if WriteVelocity
    P               = single(P_3D(:));                                                  % Pressure
    Velocity        = [single(VX_3D(:)) single(VY_3D(:))    single(VZ_3D(:))    ];      % Velocity                                       % Velocity @ every point
end

if WriteViscosity
    Eta             = single(ETA_3D(:));                                                % Viscosity
end

if WriteTCond
    TCond           = single(K_3D(:));                                                 % Thermal conductivity
end

if WriteComposition
    Composition     = single(C_3D(:));                                                % Composition
end

if WriteBasalt
    Basalt          = single(BS_3D(:));                                                % Basalt fraction
end

if WriteHPE
    HPE             = single(HPE_3D(:));                                               % Heat-producing elements
end

if WriteMeltFrac
    MeltFraction     = single(F_3D(:));                                                % Melt Fraction
end

if WriteStress
    Stress           = single(S_3D(:));                                                % Stress
end

if WriteStrainRate
    StrainRate       = single(E_3D(:));                                                % Strain rate
end

if WriteDislCreep
    RDislocation    = single(RD_3D(:));                                                % Disl. creep ratio
end
   
if WritePlasticity
    RPlasticity      = single(RP_3D(:));                                               % Plasticity fraction
end

if WritePhase
    Phasefield       = single(PH_3D(:));
end

if WriteAge
    Age = single(AGE_3D(:));
end

if WriteSmNd142
    Sm142 = single(SM142_3D(:));
    Nd142 = single(ND142_3D(:));
    Nd144 = single(ND144_3D(:));
end

if WriteHfW182
    Hf182 = single(HF182_3D(:));
    W182 = single(W182_3D(:));
    W184 = single(W184_3D(:));
end


%==========================================================================
% Definitions and initialization
sizeof_Float32  =   4;      
sizeof_Float64  =   4;     
sizeof_UInt32   =   4;     
Offset          =   0;      % Initial offset

%--------------------------------------------------------------------------
% Write the header for a structured grid:
%--------------------------------------------------------------------------
fname_vtk       = [fname,'_',num2str(1000000+fname_number),'.vts'];
fid             = fopen(fname_vtk,'w','b');           % note the 'b': not doing BigEndian does not work with MATLAB!
fprintf(fid,'<?xml version="1.0"?> \n');
fprintf(fid,'<VTKFile type="StructuredGrid" version="0.1" byte_order="BigEndian" >\n');
fprintf(fid,'  <StructuredGrid  WholeExtent="%i %i %i %i %i %i">\n', [0 size(X_3D,1)-1 0 size(X_3D,2)-1 0 size(X_3D,3)-1]);
fprintf(fid,'  <Piece Extent="%i %i %i %i %i %i">\n', [0 size(X_3D,1)-1 0 size(X_3D,2)-1 0 size(X_3D,3)-1]);

%--------------------------------------------------------------------------
% Add point-wise data
%--------------------------------------------------------------------------
fprintf(fid,'    <PointData Scalars="T" Vectors="Velocity"  >\n');

% TEMPERATURE -----------
if ASCII
    % ASCII:
    fprintf(fid,'      <DataArray type="Float32" Name="T" format="ascii">\n');
    for i=1:length(T)
        fprintf(fid,'        %g \n',single(T(i)));
    end
else
    % BINARY:
    fprintf(fid,'      <DataArray type="Float32" Name="T" format="appended" offset="%i">\n', int32(Offset));
    Offset = Offset + length(T(:))*sizeof_Float32 + 1*sizeof_UInt32;
end
fprintf(fid,'      </DataArray>\n');
% -----------------------

if WriteViscosity
    % VISCOSITY---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Viscosity" format="ascii">\n');
        for i=1:length(Eta)
            fprintf(fid,'        %g \n',single(Eta(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Viscosity" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Eta(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteDislCreep
    % R_DISLOCATION---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="R_DislCreep" format="ascii">\n');
        for i=1:length(Eta)
            fprintf(fid,'        %g \n',single(RDislocation(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="R_DislCreep" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(RDislocation(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WritePlasticity
    % R_PLASTICITY---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="R_Plasticity" format="ascii">\n');
        for i=1:length(Eta)
            fprintf(fid,'        %g \n',single(RPlasticity(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="R_Plasticity" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(RPlasticity(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end
    
if WriteTCond
    % THERMAL CONDUCTIVITY---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="TCond" format="ascii">\n');
        for i=1:length(Eta)
            fprintf(fid,'        %g \n',single(Conductivity(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="TCond" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(TCond(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteComposition
    % COMPOSITION---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Composition" format="ascii">\n');
        for i=1:length(Composition)
            fprintf(fid,'        %g \n',single(Composition(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Composition" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Composition(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteBasalt
    % BASALT FRACTION---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Basalt" format="ascii">\n');
        for i=1:length(Composition)
            fprintf(fid,'        %g \n',single(Basalt(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Basalt" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Basalt(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteHPE
    % COMPOSITION---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="HPE" format="ascii">\n');
        for i=1:length(HPE)
            fprintf(fid,'        %g \n',single(HPE(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="HPE" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(HPE(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end


if WriteMeltFrac
    % MELT FRACTION---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Melt Fraction" format="ascii">\n');
        for i=1:length(MeltFraction)
            fprintf(fid,'        %g \n',single(MeltFraction(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Melt Fraction" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(MeltFraction(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WritePhase
    % Phase Field---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Phase Field" format="ascii">\n');
        for i=1:length(Composition)
            fprintf(fid,'        %g \n',single(Phasefield(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Phase Field" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Phasefield(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end


if WriteStress
    % STRESS---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Stress" format="ascii">\n');
        for i=1:length(Stress)
            fprintf(fid,'        %g \n',single(Stress(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Stress" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Stress(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteStrainRate
    % STRAIN RATE---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Strain rate" format="ascii">\n');
        for i=1:length(StrainRate)
            fprintf(fid,'        %g \n',single(StrainRate(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Strain rate" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(StrainRate(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteMeltComp
    % MELT COMPOSITION---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Melt composition" format="ascii">\n');
        for i=1:length(MeltComposition)
            fprintf(fid,'        %g \n',single(MeltComposition(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Melt composition" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(MeltComposition(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteAge
    % Age---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Age" format="ascii">\n');
        for i=1:length(Age)
            fprintf(fid,'        %g \n',single(Age(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Age" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Age(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteSm142
    % Sm142---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Sm142" format="ascii">\n');
        for i=1:length(Sm142)
            fprintf(fid,'        %g \n',single(Sm142(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Sm142" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Sm142(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteNd142
    % Age---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Nd142" format="ascii">\n');
        for i=1:length(Nd142)
            fprintf(fid,'        %g \n',single(Nd142(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Nd142" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Nd142(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteNd144
    % Age---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Nd144" format="ascii">\n');
        for i=1:length(Nd144)
            fprintf(fid,'        %g \n',single(Nd144(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Nd144" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Nd144(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteHf182
    % Age---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Hf182" format="ascii">\n');
        for i=1:length(Hf182)
            fprintf(fid,'        %g \n',single(Hf182(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Hf182" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Hf182(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteW182
    % Age---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="W182" format="ascii">\n');
        for i=1:length(W182)
            fprintf(fid,'        %g \n',single(W182(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="W182" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(W182(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteW184
    % Age---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="W184" format="ascii">\n');
        for i=1:length(W184)
            fprintf(fid,'        %g \n',single(W184(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="W184" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(W184(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end




if WriteVelocity
    
    % PRESSURE---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="P" format="ascii">\n');
        for i=1:length(P)
            fprintf(fid,'        %g \n',single(P(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="P" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(P(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------

    % VELOCITY---------------  : VELOCITY IS A 3-component vector
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">\n');
        for i=1:length(T)
            fprintf(fid,'   %g %g %g \n',single(Velocity(i,:)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Velocity(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
    
end


fprintf(fid,'    </PointData>\n');
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

fprintf(fid,'    <Celldata>\n');
fprintf(fid,'    </Celldata>\n');


%--------------------------------------------------------------------------
% Add coordinates of structured grid 
%--------------------------------------------------------------------------
fprintf(fid,'    <Points>\n');

% ASCII
if ASCII
    fprintf(fid,'      <DataArray type="Float32" Name="Array" NumberOfComponents="3" format="ascii">\n');
    for i=1:size(Points,1)
        fprintf(fid,' %g %g %g \n',[Points(i,:)]);
    end
else
    fprintf(fid,'      <DataArray type="Float32" Name="Array" NumberOfComponents="3" format="appended" offset="%i" >\n',int32(Offset));
    
end
fprintf(fid,'      </DataArray>\n');
fprintf(fid,'    </Points>\n');
%--------------------------------------------------------------------------

fprintf(fid,'  </Piece> \n');
fprintf(fid,'  </StructuredGrid> \n');


if ~ASCII
    % Append binary data in raw format: the order in which data arrays are
    % added should be the same as how they are defined above
    fprintf(fid,'  <AppendedData encoding="raw"> \n');
    fprintf(fid,'_');

    % Add Temperature in binary format
    fwrite(fid,uint32(length(T)*sizeof_Float32),'uint32');
    fwrite(fid,single(T).'      ,   'float32');
    
    if WriteViscosity
        % Add Viscosity in binary format
        fwrite(fid,uint32(length(Eta)*sizeof_Float32),'uint32');
        fwrite(fid,single(Eta).'      ,   'float32');
    end
    
    if WriteDislCreep
        fwrite(fid,uint32(length(RDislocation)*sizeof_Float32),'uint32');
        fwrite(fid,single(RDislocation).'      ,   'float32');
    end
    
    if WritePlasticity
        fwrite(fid,uint32(length(RPlasticity)*sizeof_Float32),'uint32');
        fwrite(fid,single(RPlasticity).'      ,   'float32');
    end
    
    if WriteTCond
        fwrite(fid,uint32(length(TCond)*sizeof_Float32),'uint32');
        fwrite(fid,single(TCond).'      ,   'float32');
    end
    
    if WriteComposition
        % Add Composition in binary format
        fwrite(fid,uint32(length(Composition)*sizeof_Float32),'uint32');
        fwrite(fid,single(Composition).'      ,   'float32');
    end

    if WriteBasalt
        fwrite(fid,uint32(length(Basalt)*sizeof_Float32),'uint32');
        fwrite(fid,single(Basalt).'      ,   'float32');
    end

    if WriteHPE
        fwrite(fid,uint32(length(HPE)*sizeof_Float32),'uint32');
        fwrite(fid,single(HPE).'      ,   'float32');
    end

    if WriteMeltFrac
        % Add Melt fraction in binary format
        fwrite(fid,uint32(length(MeltFraction)*sizeof_Float32),'uint32');
        fwrite(fid,single(MeltFraction).'      ,   'float32');
    end
    
    if WritePhase
        fwrite(fid,uint32(length(Phasefield)*sizeof_Float32),'uint32');
        fwrite(fid,single(Phasefield).'      ,   'float32');
    end

    if WriteStress
        % Add Stress in binary format
        fwrite(fid,uint32(length(Stress)*sizeof_Float32),'uint32');
        fwrite(fid,single(Stress).'      ,   'float32');
    end

    if WriteStrainRate
        % Add Strain Rate in binary format
        fwrite(fid,uint32(length(StrainRate)*sizeof_Float32),'uint32');
        fwrite(fid,single(StrainRate).'      ,   'float32');
    end

    if WriteMeltComp
        % Add MeltComposition in binary format
        fwrite(fid,uint32(length(MeltComposition)*sizeof_Float32),'uint32');
        fwrite(fid,single(MeltComposition).'      ,   'float32');
    end

    if WriteAge
        % Add Age in binary format
        fwrite(fid,uint32(length(Age)*sizeof_Float32),'uint32');
        fwrite(fid,single(Age).'      ,   'float32');
    end

    if WriteNd142
        % Add Nd142 in binary format
        fwrite(fid,uint32(length(Nd142)*sizeof_Float32),'uint32');
        fwrite(fid,single(Nd142).'      ,   'float32');
    end
    if WriteNd144
        % Add Nd144 in binary format
        fwrite(fid,uint32(length(Nd144)*sizeof_Float32),'uint32');
        fwrite(fid,single(Nd144).'      ,   'float32');
    end
    if WriteSm142
        % Add Sm142 in binary format
        fwrite(fid,uint32(length(Sm142)*sizeof_Float32),'uint32');
        fwrite(fid,single(Sm142).'      ,   'float32');
    end

    if WriteHf182
        % Add Age in binary format
        fwrite(fid,uint32(length(Hf182)*sizeof_Float32),'uint32');
        fwrite(fid,single(Hf182).'      ,   'float32');
    end
    if WriteW182
        % Add W182 in binary format
        fwrite(fid,uint32(length(W182)*sizeof_Float32),'uint32');
        fwrite(fid,single(W182).'      ,   'float32');
    end
    if WriteW184
        % Add W184 in binary format
        fwrite(fid,uint32(length(W184)*sizeof_Float32),'uint32');
        fwrite(fid,single(W184).'      ,   'float32');
    end

    if WriteVelocity
        % Add Pressure in binary format
        fwrite(fid,uint32(length(P)*sizeof_Float32),'uint32');
        fwrite(fid,single(P).'      ,   'float32');

        % Add Velocity in binary format
        fwrite(fid,uint32(length(Velocity(:))*sizeof_Float32),'uint32');
        fwrite(fid,single(Velocity).' ,   'float32');
    end

    % Add Coordinates in binary format
    fwrite(fid,uint32(length(Points(:))*sizeof_Float32),'uint32');
    fwrite(fid,single(Points).' ,   'float32');

    fprintf(fid,'  </AppendedData> \n');
end


fprintf(fid,'</VTKFile>\n');
fclose(fid);

% cd(start_dir)

if ~ASCII
    disp(['Created Binary XML-VTK output file ',fname_vtk])
else
    disp(['Created ASCII XML-VTK output file ',fname_vtk])
end