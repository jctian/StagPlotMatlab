function [varargout] = ReadStag3D2(directory, fname_input, fname_number, Type)% Reads a STAG3D output file and transforms it into a MATLAB file
%
% Tested for 3D cartesian
%
%
% Syntax:
%     For scalar fields:
%
%       [X_3D, Y_3D, Z_3D, DATA_3D] = ReadStag3D(fname_input, fname_number, {'viscosity','temperature'});
%
%   For vector fields:
%
%       [X_3D, Y_3D, Z_3D, VX_3D, VY_3D, VZ_3D, P_3D] = ReadStag3D(fname_input, fname_number, {'viscosity','temperature'});
%
%
% Boris Kaus
%
% Modifications by Paul Tackley August 2010:
%  - extra fields: topography, topography self-grav, crustal
%  thickness, age, strain rate, geoid
%  - unimplemented 2Dfield stuff deleted: 3D routine works for 2D fields
%  - velocities & p now read into 4D arrays: yy split is done on output
%  - correction to the storage order of v* and p
%  - checks magic whether #components matches nval
%  - compatibility with magic=9 storing extra x & y points for vectors
%  - added missing scaling factor read for vectors
% July 2015: updated magic number stuff
% January 2021: Added basalt, HPE, conductivity, dislocation & plasticity


% start_dir = pwd;
% cd(directory);

FileFormat = 'n';           % native - default
%FileFormat = 'l';           % Little Endian
%FileFormat = 'b';           % Big    Endian

if fname_number<10000
    number_string = num2str(10000+fname_number);
    number_string(1)='0';
else
    number_string = num2str(fname_number);
end

switch Type
    case 'velocity'
        fname       = [fname_input,'_vp',number_string];
        scalardata  = false;
    case 'temperature'
        fname       = [fname_input,'_t',number_string];
        scalardata  = true;
    case 'viscosity'
        fname       = [fname_input,'_eta',number_string];
        scalardata  = true;
    case 'dislocation'
        fname       = [fname_input,'_disl',number_string];
        scalardata  = true;
    case 'plasticity'
        fname       = [fname_input,'_plast',number_string];
        scalardata  = true;
    case 'conductivity'
        fname       = [fname_input,'_k',number_string];
        scalardata  = true;
    case 'composition'
        fname       = [fname_input,'_c',number_string];
        scalardata  = true;
    case 'basalt'
        fname       = [fname_input,'_bs',number_string];
        scalardata  = true;
    case 'HPE'
        fname       = [fname_input,'_hpe',number_string];
        scalardata  = true;
    case 'melt fraction'
        fname       = [fname_input,'_f',number_string];
        scalardata  = true;
    case 'topography'
        fname       = [fname_input,'_cs',number_string];
        scalardata  = true;
    case 'topography self-grav'
        fname       = [fname_input,'_csg',number_string];
        scalardata  = true;
    case 'crustal thickness'
        fname       = [fname_input,'_cr',number_string];
        scalardata  = true;
    case 'age'
        fname       = [fname_input,'_age',number_string];
        scalardata  = true;
    case 'strain rate'
        fname       = [fname_input,'_ed',number_string];
        scalardata  = true;
    case 'stress'
        fname       = [fname_input,'_str',number_string];
        scalardata  = true;
    case 'geoid'
        fname       = [fname_input,'_g',number_string];
        scalardata  = true;
    case 'phasefield'
        fname       = [fname_input,'_ph',number_string];
        scalardata  = true;
    case 'melt composition'
        fname       = [fname_input,'_fc',number_string];
        scalardata  = true;
    otherwise
        error('Unknown property')
end
% if ~exist(fname)
%     % The file does not exist and we should stop processing data
% 
% %     cd(start_dir);
%     for i=1:10
%         varargout{i}    = -999;
%     end
% 
%     return
% end

if scalardata
    nval    =   1;      % temperature has only one value
else
    nval    =   4;      % assumed that we have a velocity-pressure file
end

fid         =   fopen(fname,'r',FileFormat);               % Open File

%==========================================================================
% READ HEADER
%==========================================================================
magic       = fread(fid,1,'int32');         % Version
if magic>8000
    INTSTR='int64'; FLTSTR='double';
    magic=magic-8000;
    fread(fid,1,'int32');
else
    INTSTR='int32'; FLTSTR='single';
end
if (magic<100 && nval>1) || (magic>300 && nval==1) % check #components
    error('wrong number of components in field')
end
magic       = mod(magic,100);
if magic>=9 && nval==4
    xyp = 1;     % extra ghost point in x & y direction
else
    xyp = 0;
end

nxtot       = fread(fid,1,INTSTR);         % nx total
nytot       = fread(fid,1,INTSTR);         % ny total
nztot       = fread(fid,1,INTSTR);         % nz total
nblocks     = fread(fid,1,INTSTR);         % # of blocks, 2 for yinyang
Aspect      = fread(fid,2,FLTSTR);        % Aspect ratio
nnx         = fread(fid,1,INTSTR);         % Number of parallel subdomains
nny         = fread(fid,1,INTSTR);         %  in the x,y,z and b directions
nnz         = fread(fid,1,INTSTR);         %
nnb         = fread(fid,1,INTSTR);         %

nz2         = nztot*2 + 1;
zg          = fread(fid,nz2,FLTSTR);      % z-coordinates

% compute nx, ny, nz and nb PER CPU
nx          =   nxtot/nnx;
ny          =   nytot/nny;
nz          =   nztot/nnz;
nb          =   nblocks/nnb;
npi         =   (nx+xyp)*(ny+xyp)*nz*nb*nval;      % the number of values per 'read' block


rcmb        = fread(fid,1,FLTSTR);
istep       = fread(fid,1,INTSTR);
time        = fread(fid,1,FLTSTR);
erupta_total= fread(fid,1,FLTSTR);
botT_val    = fread(fid,1,FLTSTR);


x           = fread(fid,nxtot,FLTSTR);      % x-coordinates
y           = fread(fid,nytot,FLTSTR);      % y-coordinates
z           = fread(fid,nztot,FLTSTR);      % z-coordinates

% read the parallel blocks
if scalardata
    DATA_3D = zeros(nxtot,nytot,nztot);
else
    scalefac= fread(fid,1,FLTSTR);            % scale factor
    VX_3D   = zeros(nxtot,nytot,nztot);         %   Vx
    VY_3D   = zeros(nxtot,nytot,nztot);         %   Vy
    VZ_3D   = zeros(nxtot,nytot,nztot);         %   Vz
    P_3D    = zeros(nxtot,nytot,nztot);         %   Pressure
end

for ibc=1:nnb       % loop over parallel subdomains
    for izc =1:nnz
        for iyc =1:nny
            for ixc =1:nnx
                data_CPU    =   fread(fid,npi,FLTSTR);      % read the data for this CPU

                % Create a 3D matrix from these data
                if scalardata
                    data_CPU_3D = reshape(data_CPU, [nx ny nz nb]);                   
                else
                    data_CPU_3D = reshape(data_CPU*scalefac, [nval nx+xyp ny+xyp nz nb]);
                end

                % Add local 3D matrix to global matrix
                if scalardata
                    % Scalar data
                    DATA_3D((ixc-1)*nx + (1:nx), (iyc-1)*ny + (1:ny), (izc-1)*nz + (1:nz),(ibc-1)*nb + (1:nb)) = data_CPU_3D;
                    
                else
                    % velocity-pressure data
                    VX_3D((ixc-1)*nx + (1:nx), (iyc-1)*ny +  (1:ny), (izc-1)*nz + (1:nz), (ibc-1)*nb + (1:nb)) = squeeze(data_CPU_3D(1,1:nx,1:ny,:,:));
                    VY_3D((ixc-1)*nx + (1:nx), (iyc-1)*ny +  (1:ny), (izc-1)*nz + (1:nz), (ibc-1)*nb + (1:nb)) = squeeze(data_CPU_3D(2,1:nx,1:ny,:,:));
                    VZ_3D((ixc-1)*nx + (1:nx), (iyc-1)*ny +  (1:ny), (izc-1)*nz + (1:nz), (ibc-1)*nb + (1:nb)) = squeeze(data_CPU_3D(3,1:nx,1:ny,:,:));
                    P_3D( (ixc-1)*nx + (1:nx), (iyc-1)*ny +  (1:ny), (izc-1)*nz + (1:nz), (ibc-1)*nb + (1:nb)) = squeeze(data_CPU_3D(4,1:nx,1:ny,:,:));
                end


            end
        end
    end
end

fclose(fid);                                % close file



[Y_3D, X_3D, Z_3D]  = meshgrid(y,x,z);

% prepare output data
varargout{1} = nblocks;
varargout{2} = X_3D;
varargout{3} = Y_3D;
varargout{4} = Z_3D;

if nblocks==1
    % no ying-yang
    switch Type
        case 'velocity'
            varargout{5} = VX_3D;
            varargout{6} = VY_3D;
            varargout{7} = VZ_3D;
            varargout{8} = P_3D;
            varargout{9} = time;
        otherwise
            varargout{5} = DATA_3D;
            varargout{6} = time;
            varargout{7} = rcmb;
    end
else
    % ying-yang grid
    switch Type
        case 'velocity'
            varargout{5} = VX_3D(:,:,:,1);
            varargout{6} = VY_3D(:,:,:,1);
            varargout{7} = VZ_3D(:,:,:,1);
            varargout{8} = P_3D (:,:,:,1);

            varargout{9}  = VX_3D(:,:,:,2);
            varargout{10} = VY_3D(:,:,:,2);
            varargout{11} = VZ_3D(:,:,:,2);
            varargout{12} = P_3D (:,:,:,2);

            varargout{13} = time;
        otherwise
            varargout{5} = DATA_3D(:,:,:,1);
            varargout{6} = DATA_3D(:,:,:,2);
            varargout{7} = time;
            varargout{8} = rcmb;
    end

end



% cd(start_dir);



