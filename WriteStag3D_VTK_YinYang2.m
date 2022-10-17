% Write STAG_3D format from matlab in XML-VTK format, using ASCII or BINARY file
% format, in case we are dealing with a YING-YANG GRID.
%
%
%
% Boris Kaus, 26 Feb. 2008
% 8.13.2008  Bugfix that caused compositional fields etc. not to be written
%            correctly to file

% Paul Tackley 3.08.2010: corrections to pressure and velocity writing
%                       modified corner deletion for new coordinates

% The ying-yang grid consists of two blocks, with overlapping corners.
% The corners are removed, and the grid is triangulated. A single unstructured
% mesh is generated from this, that has VT_WEDGE elements. 

ASCII           = logical(0);                       % ASCII or BINARY?

% Change to correct directory
start_dir       = pwd;
cd(directory);

%==========================================================================
% 1) Take the surface of the 2 grids, patch together and triangulate
%==========================================================================
xs1                 =   X_3D_1(:,:,end);
ys1                 =   Y_3D_1(:,:,end);
zs1                 =   Z_3D_1(:,:,end);

xs2                 =   X_3D_2(:,:,end);
ys2                 =   Y_3D_2(:,:,end);
zs2                 =   Z_3D_2(:,:,end);

r1                  =   sqrt(xs1.^2 + ys1.^2 + zs1.^2);
theta1              =   atan2(sqrt(xs1.^2 + ys1.^2),zs1);
phi1                =   atan2(ys1,xs1);


r2                  =   sqrt(xs2.^2 + ys2.^2 + zs2.^2);
theta2              =   atan2(sqrt(xs2.^2 + ys2.^2),zs2);
phi2                =   atan2(ys2,xs2);


% Cut off the corners from grid #1, which seems to do #2 as well (PJT):
theta12             = acos(sin(theta1).*sin(phi1));  % theta coords of grid 1 in coord system of grid 2
ind_corner          = find( (theta12>pi/4 & phi1>pi/2) | (theta12<3*pi/4 & phi1<-pi/2 ) );

% Form indices of remaining 
ind                           = 1:prod(size(phi1));
ind([ind_corner])=[];

% Create an array with unique r, theta & phi values
numberYin = ones(size(r1));
numberYin(find(r1==r1))           =   find(r1==r1);
numberYang = ones(size(r2));
numberYang(find(r2==r2))          =   find(r2==r2)+ max(numberYin(:));

%Closed Yin grid, w/o corners and completely closed
r                   =   [r1(ind),      ]';
theta               =   [theta1(ind)   ]';
phi                 =   [phi1(ind),    ]';
NumYin              =   [numberYin(ind)]';

% x,y,z coordinates of complete grid:
x                   =   [xs1(ind) -xs1(ind)];
y                   =   [ys1(ind)  zs1(ind)];
z                   =   [zs1(ind)  ys1(ind)];
NumYang             =   [NumYin + max(numberYin(:)) ];
Number              =   [NumYin(:);    NumYang(:)  ]';

%   remove redundant coordinates
tri                 =   convhulln([x(:), y(:),z(:)]);       % simple way to grid it
triYingYang         =   Number(tri);

x                   =   [xs1(:); xs2(:)];
y                   =   [ys1(:); ys2(:)];
z                   =   [zs1(:); zs2(:)];
%==========================================================================
% triYingYang now contains the numbers of all triangles 
%==========================================================================

%==========================================================================
% 2) Create a 3D grid with tetrahedron elements
%==========================================================================

% Number all gridpoints we have
NUMBER_1                 = ones(size(X_3D_1));
NUMBER_2                 = ones(size(X_3D_2));
NUMBER_1(find(NUMBER_1)) = find(NUMBER_1);
NUMBER_2(find(NUMBER_2)) = find(NUMBER_2) + max(NUMBER_1(:));


% Make a loop over all levels
ElementNumbers      = [];
for iz=1:size(X_3D_2,3)-1

    num_upper1   	= NUMBER_1(:,:,iz+1);
    num_upper2      = NUMBER_2(:,:,iz+1);
    num_upper       = [num_upper1(:); num_upper2(:)];

    num_lower1      = NUMBER_1(:,:,iz);
    num_lower2      = NUMBER_2(:,:,iz);
    num_lower       = [num_lower1(:); num_lower2(:)];
    ElementNumbers  = [ElementNumbers; num_upper(triYingYang), num_lower(triYingYang)];
end


%--------------------------------------------------------------------------
% Convert data into correct vector format
%--------------------------------------------------------------------------
Points                  = zeros(max(NUMBER_2(:)),3);
Points(NUMBER_1(:),1)       = X_3D_1(:);
Points(NUMBER_2(:),1)       = X_3D_2(:);
Points(NUMBER_1(:),2)       = Y_3D_1(:);
Points(NUMBER_2(:),2)       = Y_3D_2(:);
Points(NUMBER_1(:),3)       = Z_3D_1(:);
Points(NUMBER_2(:),3)       = Z_3D_2(:);

Radius                      = sqrt(sum(Points.^2,2));


Temperature              = zeros(max(NUMBER_2(:)),1);
Temperature(NUMBER_1(:)) = T_3D_1(:);
Temperature(NUMBER_2(:)) = T_3D_2(:);

if WriteVelocity
    Velocity_vec                = zeros(max(NUMBER_2(:)),3);
    Velocity_vec(NUMBER_1(:),1) = VX_3D_1(:);
    Velocity_vec(NUMBER_2(:),1) = VX_3D_2(:);
    Velocity_vec(NUMBER_1(:),2) = VY_3D_1(:);
    Velocity_vec(NUMBER_2(:),2) = VY_3D_2(:);
    Velocity_vec(NUMBER_1(:),3) = VZ_3D_1(:);
    Velocity_vec(NUMBER_2(:),3) = VZ_3D_2(:);
    
    Pressure_vec    = zeros(max(NUMBER_2(:)),1);
    Pressure_vec(NUMBER_1(:),1) = P_3D_1(:);
    Pressure_vec(NUMBER_2(:),1) = P_3D_2(:);
end

if WriteViscosity
    Viscosity_vec                = zeros(max(NUMBER_2(:)),1);
    Viscosity_vec(NUMBER_1(:),1) = ETA_3D_1(:);
    Viscosity_vec(NUMBER_2(:),1) = ETA_3D_2(:);
end

if WriteDislCreep
    DislCreep_vec                = zeros(max(NUMBER_2(:)),1);
    DislCreep_vec(NUMBER_1(:),1) = RD_3D_1(:);
    DislCreep_vec(NUMBER_2(:),1) = RD_3D_2(:);
end

if WritePlasticity
    Plasticity_vec                = zeros(max(NUMBER_2(:)),1);
    Plasticity_vec(NUMBER_1(:),1) = RP_3D_1(:);
    Plasticity_vec(NUMBER_2(:),1) = RP_3D_2(:);
end

if WriteTCond
    TCond_vec                = zeros(max(NUMBER_2(:)),1);
    TCond_vec(NUMBER_1(:),1) = K_3D_1(:);
    TCond_vec(NUMBER_2(:),1) = K_3D_2(:);
end

if WriteComposition
    Composition_vec                = zeros(max(NUMBER_2(:)),1);
    Composition_vec(NUMBER_1(:),1) = C_3D_1(:);
    Composition_vec(NUMBER_2(:),1) = C_3D_2(:);
end

if WriteBasalt
    Basalt_vec                = zeros(max(NUMBER_2(:)),1);
    Basalt_vec(NUMBER_1(:),1) = BS_3D_1(:);
    Basalt_vec(NUMBER_2(:),1) = BS_3D_2(:);
end

if WriteHPE
    HPE_vec                = zeros(max(NUMBER_2(:)),1);
    HPE_vec(NUMBER_1(:),1) = HPE_3D_1(:);
    HPE_vec(NUMBER_2(:),1) = HPE_3D_2(:);
end

if WriteMeltFrac
    MeltFraction_vec                = zeros(max(NUMBER_2(:)),1);
    MeltFraction_vec(NUMBER_1(:),1) = F_3D_1(:);
    MeltFraction_vec(NUMBER_2(:),1) = F_3D_2(:);
end

if WriteStress
    Stress_vec                = zeros(max(NUMBER_2(:)),1);
    Stress_vec(NUMBER_1(:),1) = S_3D_1(:);
    Stress_vec(NUMBER_2(:),1) = S_3D_2(:);
end

if WriteStrainRate
    StrainRate_vec                = zeros(max(NUMBER_2(:)),1);
    StrainRate_vec(NUMBER_1(:),1) = E_3D_1(:);
    StrainRate_vec(NUMBER_2(:),1) = E_3D_2(:);
end


%==========================================================================
% 3) Write VTK file (unstructured mesh)
%==========================================================================
ElementNumbers  = ElementNumbers-1;      % VTK is zero-based 


%==========================================================================
% Definitions and initialization
sizeof_Float32  =   4;      
sizeof_Float64  =   4;     
sizeof_UInt32   =   4; 
sizeof_UInt8    =   1; 

Offset          =   0;      % Initial offset

%--------------------------------------------------------------------------
% Write the header for a structured grid:
%--------------------------------------------------------------------------
fname_vtk       = [fname,'_',num2str(1000000+fname_number),'.vtu'];
fid             = fopen(fname_vtk,'w','b');           % note the 'b': not doing BigEndian does not work with MATLAB!
fprintf(fid,'<?xml version="1.0"?> \n');
fprintf(fid,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian" >\n');
fprintf(fid,'  <UnstructuredGrid>\n');
fprintf(fid,'    <Piece NumberOfPoints="%i"  NumberOfCells="%i">\n', int64(size(Points,1)), int64(size(ElementNumbers,1)));



%--------------------------------------------------------------------------
% Add point-wise data
%--------------------------------------------------------------------------
fprintf(fid,'    <PointData Scalars="T" Vectors="Velocity"  >\n');

% TEMPERATURE -----------
if ASCII
    % ASCII:
    fprintf(fid,'      <DataArray type="Float32" Name="T" format="ascii">\n');
    for i=1:length(Temperature)
        fprintf(fid,'        %g \n',single(Temperature(i)));
    end
else
    % BINARY:
    fprintf(fid,'      <DataArray type="Float32" Name="T" format="appended" offset="%i">\n', int64(Offset));
    Offset = Offset + length(Temperature(:))*sizeof_Float32 + 1*sizeof_UInt32;
end
fprintf(fid,'      </DataArray>\n');
% -----------------------


if WriteViscosity
    % VISCOSITY---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Viscosity" format="ascii">\n');
        for i=1:length(Viscosity_vec)
            fprintf(fid,'        %g \n',single(Viscosity_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Viscosity" format="appended" offset="%i">\n',int64(Offset));
        Offset = Offset + length(Viscosity_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteDislCreep
    % DISLOCATION CREEP RATIO---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="R_DislCreep" format="ascii">\n');
        for i=1:length(DislCreep_vec)
            fprintf(fid,'        %g \n',single(DislCreep_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="R_DislCreep" format="appended" offset="%i">\n',int64(Offset));
        Offset = Offset + length(R_DislCreep_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WritePlasticity
    % PLASTICITY FRACTION---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="R_Plasticity" format="ascii">\n');
        for i=1:length(Plasticity_vec)
            fprintf(fid,'        %g \n',single(Plasticity_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Plasticity" format="appended" offset="%i">\n',int64(Offset));
        Offset = Offset + length(Plasticity_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteTCond
    % THERMAL CONDUCTIVITY---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="TCond" format="ascii">\n');
        for i=1:length(TCond_vec)
            fprintf(fid,'        %g \n',single(TCond_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="TCond" format="appended" offset="%i">\n',int64(Offset));
        Offset = Offset + length(TCond_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteComposition
    % COMPOSITION---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Composition" format="ascii">\n');
        for i=1:length(Composition_vec)
            fprintf(fid,'        %g \n',single(Composition_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Composition" format="appended" offset="%i">\n',int64(Offset));
        Offset = Offset + length(Composition_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteBasalt
    % BASALT FRACTION---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Basalt" format="ascii">\n');
        for i=1:length(Basalt_vec)
            fprintf(fid,'        %g \n',single(Basalt_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Basalt" format="appended" offset="%i">\n',int64(Offset));
        Offset = Offset + length(Basalt_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteHPE
    % HEAT-PRODUCING ELEMENTS---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="HPE" format="ascii">\n');
        for i=1:length(HPE_vec)
            fprintf(fid,'        %g \n',single(HPE_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="HPE" format="appended" offset="%i">\n',int64(Offset));
        Offset = Offset + length(HPE_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteMeltFrac
    % MELT FRACTION---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Melt fraction" format="ascii">\n');
        for i=1:length(MeltFraction_vec)
            fprintf(fid,'        %g \n',single(MeltFraction_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Melt fraction" format="appended" offset="%i">\n',int64(Offset));
        Offset = Offset + length(MeltFraction_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteStress
    % STRESS---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Stress" format="ascii">\n');
        for i=1:length(Stress_vec)
            fprintf(fid,'        %g \n',single(Stress_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Stress" format="appended" offset="%i">\n',int64(Offset));
        Offset = Offset + length(Stress_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteStrainRate
    % STRAIN RATE---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Strain rate" format="ascii">\n');
        for i=1:length(StrainRate_vec)
            fprintf(fid,'        %g \n',single(StrainRate_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Strain rate" format="appended" offset="%i">\n',int64(Offset));
        Offset = Offset + length(StrainRate_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteVelocity
    
    % PRESSURE---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="P" format="ascii">\n');
        for i=1:length(Pressure_vec)
            fprintf(fid,'        %g \n',single(Pressure_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="P" format="appended" offset="%i">\n',int64(Offset));
        Offset = Offset + length(Pressure_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------

    % VELOCITY---------------  : VELOCITY IS A 3-component vector
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">\n');
        for i=1:length(T)
            fprintf(fid,'   %g %g %g \n',single(Velocity_vec(i,:)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="appended" offset="%i">\n',int64(Offset));
        Offset = Offset + length(Velocity_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
    
end
fprintf(fid,'    </PointData>\n');

%--------------------------------------------------------------------------
% Add coordinates of structured grid 
%--------------------------------------------------------------------------
fprintf(fid,'    <Points>\n');

% ASCII
if ASCII
    fprintf(fid,'      <DataArray type="Float32" Name="Array" NumberOfComponents="3" format="ascii">\n');
    for i=1:size(Points,1)
        fprintf(fid,'         %g %g %g \n',[Points(i,:)]);
    end
else
     fprintf(fid,'      <DataArray type="Float32" Name="Array" NumberOfComponents="3" format="appended" offset="%i" >\n',int64(Offset));
     Offset = Offset + length(Points(:))*sizeof_Float32 + 1*sizeof_UInt32;
end
fprintf(fid,'      </DataArray>\n');
fprintf(fid,'    </Points>\n');


%--------------------------------------------------------------------------
% Add CELLS data
%--------------------------------------------------------------------------
fprintf(fid,'    <Cells>\n');

% Connectivity -----------
if ASCII
    % ASCII:
    fprintf(fid,'      <DataArray type="Int32" Name="connectivity" format="ascii">\n');
    for i=1:size(ElementNumbers,1)
        fprintf(fid,'        %i %i %i %i %i %i \n',int64(ElementNumbers(i,:)));
    end
    
else
   fprintf(fid,'      <DataArray type="Int32" Name="connectivity" format="appended" offset="%i">\n',int64(Offset));
   Offset = Offset + length(ElementNumbers(:))*sizeof_UInt32 + 1*sizeof_UInt32;
    
end
fprintf(fid,'      </DataArray>\n');

% Offsets -----------
offsets = cumsum(ones(size(ElementNumbers,1),1)*6);
if ASCII
    % ASCII:
    fprintf(fid,'  <DataArray type="Int32" Name="offsets" format="ascii">\n');
    for i=1:size(ElementNumbers,1)
        fprintf(fid,'        %i\n',int64(offsets(i)));
    end

else
    fprintf(fid,'      <DataArray type="Int32" Name="offsets" format="appended" offset="%i">\n',int64(Offset));
    Offset = Offset + length(offsets(:))*sizeof_UInt32 + 1*sizeof_UInt32;
end
fprintf(fid,'      </DataArray>\n');

% types -----------
types = ones(size(ElementNumbers,1),1)*13;
if ASCII
    % ASCII:
    fprintf(fid,'      <DataArray type="UInt8" Name="types" format="ascii">\n');
    for i=1:size(ElementNumbers,1)
        fprintf(fid,'        %i\n',uint8(13));
    end

else
    fprintf(fid,'      <DataArray type="UInt8" Name="types" format="appended" offset="%i">\n',int64(Offset));
    Offset = Offset + length(types(:))*sizeof_UInt8 + 1*sizeof_UInt32;
end
fprintf(fid,'      </DataArray>\n');


fprintf(fid,'    </Cells>\n');
% -----------------------


fprintf(fid,'    </Piece>\n');
fprintf(fid,'  </UnstructuredGrid>\n');


%--------------------------------------------------------------------------

if ~ASCII
    % Append binary data in raw format: the order in which data arrays are
    % added should be the same as how they are defined above
    fprintf(fid,'  <AppendedData encoding="raw"> \n');
    fprintf(fid,'_');

    % Add Temperature in binary format
    fwrite(fid,uint32(length(Temperature(:))*sizeof_Float32),'uint32');
    fwrite(fid,single(Temperature(:)).'          ,   'float32');
    
    if WriteViscosity
        % Add Viscosity in binary format
        fwrite(fid,uint32(length(Viscosity_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Viscosity_vec(:)).'     ,   'float32');
    end

    if WriteDislCreep
        fwrite(fid,uint32(length(DislCreep_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(DislCreep_vec(:)).'     ,   'float32');
    end

    if WritePlasticity
        fwrite(fid,uint32(length(Plasticity_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Plasticity_vec(:)).'     ,   'float32');
    end

    if WriteTCond
        fwrite(fid,uint32(length(TCond_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(TCond_vec(:)).'     ,   'float32');
    end

    if WriteComposition
        % Add Composition in binary format
        fwrite(fid,uint32(length(Composition_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Composition_vec(:)).'   ,   'float32');
    end
    
    if WriteBasalt
        fwrite(fid,uint32(length(Basalt_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Basalt_vec(:)).'   ,   'float32');
    end
    
    if WriteHPE
        fwrite(fid,uint32(length(HPE_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(HPE_vec(:)).'   ,   'float32');
    end
    
    if WriteMeltFrac
        % Add Melt fraction in binary format
        fwrite(fid,uint32(length(MeltFraction_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(MeltFraction_vec(:)).'  ,   'float32');
    end

    if WriteStress
        % Add Stress in binary format
        fwrite(fid,uint32(length(Stress_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Stress_vec(:)).'  ,   'float32');
        end

    if WriteStrainRate
        % Add Strain rate in binary format
        fwrite(fid,uint32(length(StrainRate_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(StrainRate_vec(:)).'  ,   'float32');
    end

    if WriteVelocity
        % Add Pressure in binary format
        fwrite(fid,uint32(length(Pressure_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Pressure_vec(:)).'      ,   'float32');

        % Add Velocity in binary format
        fwrite(fid,uint32(length(Velocity_vec(:))*sizeof_Float32),'uint32');
        fwrite(fid,single(Velocity_vec).' ,   'float32');
    end

    % Add Coordinates in binary format
    fwrite(fid,uint32(length(Points(:))*sizeof_Float32),'uint32');
    fwrite(fid,single(Points).' ,   'float32');
    
    % Add Element connectivity in binary format
    fwrite(fid,uint32(length(ElementNumbers(:))*sizeof_UInt32),'uint32');
    fwrite(fid,uint32(ElementNumbers).' ,   'int32');

    % Add offsets array in binary format
    fwrite(fid,uint32(length(offsets(:))*sizeof_UInt32),'uint32');
    fwrite(fid,uint32(offsets).' ,   'int32');
    
    % Add types array in binary format
    fwrite(fid,uint32(length(types(:))*sizeof_UInt32),'uint32');
    fwrite(fid,uint8(types).' ,   'int8');
    
    fprintf(fid,'  </AppendedData> \n');
end
fprintf(fid,'</VTKFile>\n');

fclose(fid);

if ~ASCII
    disp(['Created Binary YinYang XML-VTK output file ',fname_vtk])
else
    disp(['Created ASCII  YinYang XML-VTK output file ',fname_vtk])
end

cd(start_dir);                                  % Change back to original directory

