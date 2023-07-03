function [data, header] = readrec_V4_2(filename, varargin)
% [data, header] = readrec_V4_2(filename[, 'quiet'][, 'par'])
%
%  This will read in a par/rec file from the Philips scanner.
%  It will read in REC files even with real and imaginary data
%  as well as do the appropriate signal intensity adjustments
%  as the REC file stores uint16, but the underlying data can be
%  positive or negative.
%
% If filename is left off call, GUI file chooser is used.
% If filename is a struct, it assumes this is a par struct returned from
%    a previous call and will read only the rec.
% If 'quiet' is provided, operates silently without displaying header
%  info or progress bar.
% If 'par' is provided, only par file is read and only header struct is
%  returned.
%

% Craig Jones (craig@mri.jhu.edu)  April 28, 2004
% 20040510 CJ - changed reading loop to be over size(A,1)
% 20040527 CJ - checked for 'rec' and 'REC'
% 20040624 CJ - fixed TeX interpreter,
% 20040708 CJ - correctly reads in types for 1.5T data
% 20040711 CJ - fixed number of dynamics
% 20040916 CJ - added in reading of echoes
% 20050503 CJ - fixed for reading version 3/4 PAR files.
% 20050520 CJ - check for existence of par file.
% 20050606 CJ - added for cardiac phases
% 20060303 CJ - Fixed problem with V in patient name
% 20060308 CJ - Fixed Jonathan's problems with passing in a par file
% Mina Kim (mina@mri.jhu.edu)    March 26,  2007
% 20070326 MK, JSG - Fixed for reading version V4.1 par file
% Alan Huang (Alan.J.Huang@gmail.com)   November, 2007
% 20071105 AH, JH - Fixed for reading version V4.2 par file
% 20071106 AH, JF, IL - Fixed for reading ASL images
% 20120327 JSG - add uigetfile when no filename
% 20130326 JSG - add par only read and rec only read, expand header
%    structure to include all fields from par file
% 20150814 PG - added option to skip voxel intensity scaling (code by Wouter/Eva)
%

%% args processing
data = []; header = [];
noscale = 0;
quiet = 0;
par = 0;
i = 1;
while i <= length(varargin)
    switch varargin{i}
        case 'noscale'
            noscale = 1;
        case 'quiet'
            quiet = 1;
        case 'par'
            par = 1;
        otherwise
            error ('Unrecognized property : %s', varargin{i});
    end
    i = i + 1;
end

%% filename processing
if (nargin > 0) && isstruct(filename)
    % par header struct passed in
    header = filename;
    filename_par = header.filename;
    [PathName, FileName, ~] = fileparts(filename_par);
    filename = fullfile (PathName, [FileName '.rec']);
else
    % use GUI to get filename
    if (nargin == 0) || (strcmp (filename, ''))
        [FileName,PathName] = uigetfile ('*.par;*.PAR', 'Select a PAR file');
        if isequal(FileName,0)
            return;
        end
        filename = fullfile(PathName, FileName);
    end
    % split up fileame - try to get case correct
    [PathName, FileName, FileExt] = fileparts(filename);
    if ( strcmp (FileExt, '.par' ) )
        filename_par = filename;
        filename = fullfile (PathName, [FileName '.rec']);
    elseif ( strcmp(FileExt, '.PAR' ) )
        filename_par = filename;
        filename = fullfile (PathName, [FileName '.REC']);
    elseif ( strcmp(FileExt, '.rec') )
        filename_par = fullfile (PathName, [FileName '.par']);
    elseif ( strcmp(FileExt, '.REC') )
        filename_par = fullfile (PathName, [FileName '.PAR']);
    else
        error('Yo dude, %s does not seem to be a PAR or REC file', filename);
    end
    
    % load in the PAR file
    header = parseHeader(filename_par, quiet);
    if header.nrows == 0
        return;
    end
end

%% if par is true, return header only - no rec file read
if par == 1
    return
end

%%
A = header.tbl;
%if (header.par_version >= 4)
%    if (header.corbadump)
%        typeind = 5;
%    else
%        typeind = 6; %??need to check this
%    end
%else
typeind = 5;
%end

% Number of types of scans:  0 = Magnitude, 1 = Real, 2 = Imaginary,
%   16 = corbadump
types =    unique(A(:,typeind));
dynamics = unique(A(:,3));
echoes =   unique(A(:,2));
phases =   unique(A(:,4));

if (header.par_version > 4)
    bdir = unique(A(:,43));
    if (header.par_version == 4.2)
        ASL_type = unique(A(:,49));
    end
end
iri = header.tblcols.rescale_int;
irs = header.tblcols.rescale_slope;
iss = header.tblcols.scale_slope;


%%  open rec file
if (strfind(filename, '.gz') )
    tmpname = tempname;
    unix(sprintf('gzcat %s > %s', filename, tmpname));
else
    tmpname = filename;
end

fp = fopen(tmpname, 'rb', 'l');
if (fp == -1 )
    error('readrec: file %s does not exist', tmpname);
end

% progress bar
if ~quiet
    h = waitbar(0, strrep(sprintf('Reading in %s', filename),'_','_'));
end

%% read rec file
s = size(A,1);
sz = header.sizes;
sz(1) = header.sizes(2); % rows and cols are swapped when read
sz(2) = header.sizes(1);
if noscale
    data = zeros(sz, 'uint16');
else
    data = zeros(sz);
end

% eric, non-square image flag
isNonSquare = sz(1) ~= sz(2);

% when file is sorted (as is usual) ii = idx+1
% indexes go from 0 to n-1, n = # of images in file
for idx = 0:s-1
    % use index column to find info for next data block to read
    ii = find(A(:,7)==idx);
    
    %  Determine the rescale slope, intercept and other scaling factors.
    %     intercept = A(ii,iri) / A(ii,irs); % ri is usually 0 for magn data
    %     slope = 1.0 / A(ii,iss);
    
    % for Qflow
    slope = A(ii,irs);
    intercept = A(ii,iri);
    
    if noscale
        d = fread (fp, [header.nrows, header.ncols], 'uint16');
    else
        d = fread (fp, [header.nrows, header.ncols], 'int16') * slope + intercept;
    end
    if isNonSquare
        d = d';
    end
    
    if (size(d) == header.sizes(2:-1:1))
        if (header.par_version == 4.2)
            % v 4.2 - diffusion & ASL type dimension
            data(:,:, A(ii,1), echoes==A(ii,2), dynamics==A(ii,3),...
                types==A(ii,typeind), phases==A(ii,4), ...
                bdir==A(ii,43), ASL_type==A(ii,49)) = d;
        elseif (header.par_version >= 4)
            % v 4.1 - no ASL type dimension
            data(:,:, A(ii,1), find(echoes==A(ii,2)), find(dynamics==A(ii,3)),...
                find(types==A(ii,typeind)), find(phases==A(ii,4)), ...
                find(bdir==A(ii,43))) = d; %#ok<FNDSB>
        else
            % v3 - no ASL, diffusion shows up as dynamic
            data(:,:, A(ii,1), find(echoes==A(ii,2)), find(dynamics==A(ii,3)),...
                find(types==A(ii,typeind)), find(phases==A(ii,4))) = d; %#ok<FNDSB>
        end
    else
        if ~quiet
            error ('data size error');
        end
    end
    if ~quiet
        waitbar (idx/s, h)
    end
end

if isNonSquare
    data = permute(data,[2 1 3 4 5 6 7 8 9]);
    temp = header.nrows;
    header.nrows = header.ncols;
    header.ncols = temp;
end
fclose(fp);

%% reorder rows and cols
if ~quiet
    waitbar(1, h, 'Re-Ordering...');
end
data = permute(data, [2 1 3:length(size(data))]);

if( strfind(filename, '.gz') )
    delete(tmpname);
end
if ~quiet
    close(h);
end

%========================================================================
function header = parseHeader(filename_par, quiet)
%  parseHeader - Parse a Philips PAR file for some important parameters
%

%  Craig Jones (craig@mri.jhu.edu)
%  20040616 - fixed up reading in the header for stupid Windoze boxes

%  Alan Huang (Alan.J.Huang@gmail.com)
%  20071106 - added nASL_type parameter for patch 2.5 (V4.2)

% jsg - add all par header params to the struct

header.nrows = 0; header.ncols = 0; header.nslices = 0; header.nechoes = 0;
header.ndynamics = 0; header.nphases = 0; header.par_version = 0;
header.corbadump = 0;
header.filename = filename_par;

fp = fopen(filename_par, 'rt');

if( fp == -1 )
    fprintf(1,'readrec: file %s does not exist\n',  filename_par);
    return;
end

%% first section
fgetl(fp); %# === DATA DESCRIPTION FILE ====
fgetl(fp); %#
fgetl(fp); %# CAUTION - Inves...
fgetl(fp); %# Limited by Federal Law ...
fgetl(fp); %#
line = fgetl(fp); %# Dataset name: <filename>
if isa(line, 'char')
    s = regexpi(line, 'Dataset name: (?<data>.*)', 'names');
    if (isempty(s))
        fprintf (1,'readrec: Bad PAR file1: %s\n', filename_par);
        return
    end
else
    fprintf (1,'readrec: Bad PAR file2: %s\n', filename_par);
    return
end
header.dataset = s(1).data;
% flag for par made by real-time CorbaDumper program
header.corbadump = strncmp ('Dump-', header.dataset, 1) == 1;
fgetl(fp); %#
if isa(line, 'char')
    line = fgetl(fp); %# CLINICAL TRYOUT  Research image export tool   V4.2
    s = regexpi(line, 'V(?<ver>[0-9](\.[0-9]){0,1})', 'names');
    if (isempty(s))
        fprintf (1,'readrec: Bad PAR file3 %s\n', filename_par);
        return
    end
else
    fprintf (1,'readrec: Bad PAR file4 %s\n', filename_par);
    return
end
header.par_version = str2double(s(1).ver);
% fail if par version not recognized
if ((header.par_version ~= 3) && (header.par_version ~= 4) &&...
        (header.par_version ~= 4.1) && (header.par_version ~= 4.2))
    fprintf (1,'readrec: PAR file version not recognized :%s\n',...
        header.par_version);
    return
end
fgetl(fp); %#
fgetl(fp); %# === GENERAL INFORMATION ----
fgetl(fp); %#

%% second section - lines begin .
pars = {
    %label                struct id    isnum
    'Patient name',       'patname',     0;
    'Examination name',   'examname',    0;
    'Protocol name',      'protocol',    0;
    'Examination date',   'datetime',    0;
    'Series Type',        'seriestype',  0; %v4
    'Acquisition nr',     'acq',         1;
    'Reconstruction nr',  'recon',       1;
    'Scan Duration',      'duration',    1;
    'number of cardiac',  'nphases',     1;
    'number of echoes',   'nechoes',     1;
    'number of slices',   'nslices',     1;
    'number of dynamics', 'ndynamics',   1;
    'number of mixes',    'nmixes',      1;
    'Image pixel size',   'pixsize',     1; %v3
    'Patient position',   'position',    0; %v4
    'Preparation direct', 'prepdir',     0; %v4
    'Technique',          'technique',   0;
    'Scan resolution',    'scanres',     1;
    'Scan mode',          'scanmode',    0;
    'Scan percentage',    'scanpct',     1; %v3
    'Recon resolution',   'reconres',    1; %v3
    'Number of aver',     'nav',         1; %v3
    'Repetition time',    'tr',          1;
    'FOV (',              'fov',         1;
    'Slice thickness',    'slthick',     1; %v3
    'Slice gap',          'slgap',       1; %v3
    'Water Fat shift',    'wfs',         1;
    'Angulation ',        'angulation',  1;
    'Off Centre ',        'off_center',  1;
    'Flow compensation',  'flowcomp',    1;
    'Presaturation',      'presat',      1;
    'Cardiac freq',       'cardfreq',    1; %v3
    'Min. RR',            'minrr',       1; %v3
    'Max. RR',            'maxrr',       1; %v3
    'Phase encoding vel', 'pevelocity',  1;
    'MTC',                'mtc',         1;
    'SPIR',               'spir',        1;
    'EPI factor',         'epifactor',   1;
    'TURBO factor',       'turbo',       1; %v3
    'Dynamic scan',       'dynamic',     1;
    'Diffusion  ',        'diffusion',   1;
    'Diffusion echo',     'diffecho',    1;
    'diffusion values',   'diffvalues',  1; %v4
    'gradient orients',   'gradorients', 1; %v4
    'Inversion delay',    'ti',          1; %v3
    'label types',        'nASL_type',   1  %v4.2
    };

start = 1;
while( 1 )
    line = fgetl(fp);
    % finished if line begins # sl ec - heading of the 3rd section
    if (regexpi(line, '^#\s*sl\s+ec'))
        break
    end
    
    for i = start:size(pars,1)
        % if (getPar (line, pars{i,1}, pars{i,2}, pars{i,3}))
        if (strfind(upper(line), upper(pars{i,1})) > 0)  % paul made this case-insensitive
            value = line(strfind(line, ':')+1:end);
            if (pars{i,3} == 1)
                value = str2num(value); %#ok<ST2NM>
            else
                value = strtrim(value);
            end
            header.(pars{i,2}) = value;
            start = i + 1;
            break;
        end
    end
    
end

fgetl(fp); % skip blank line following heading
header.annotation = 'ap fh rl'; % ???

%  Look for number of rows and columns for v3
if (header.par_version < 4)
    header.nrows = header.reconres(1);
    header.ncols = header.reconres(2);
end

%% third section - table values
A = [];
ii = 1;
while( 1 )
    line = fgetl(fp);
    % empty line or comment at end
    if ((length(line) < 2) || (strncmp('# =', line, 3) == 1))
        break;
    end
    line = strtrim(line);
    A(ii,:) = str2num(line); %#ok<AGROW,ST2NM>
    ii = ii + 1;
end
fclose(fp);

% set these based on table rows rather than second header
header.ndynamics = length(unique(A(:,3)));
header.nechoes = length(unique(A(:,2)));
if (header.par_version == 4.2)
    header.nASL_type = length(unique(A(:,49)));
else
    header.nASL_type = [];
end
%  Added for the V4+ PAR files - assume same for every row
if (header.par_version >= 4 )
    header.nrows = A(1,10);
    header.ncols = A(1,11);
    
             header.pixdim = [A(1,29) A(1,30) A(1,23)+A(1,24)];
    % for PROUD, overwrite to true extent (FOV) / matrix size
%     header.pixdim = [header.fov(2)/header.ncols header.fov(1)/header.nrows header.fov(3)/header.nslices];
    
end

%% additional header fields
% copy bottom table into header
header.tbl = A;

% labels for tbl columns
header.tblcols = {'slice' 'echo' 'dynamic' 'card_phase' 'img_type' 'seq' ...
    'index'};
if (header.par_version > 4)
    header.tblcols = [header.tblcols {'pix_size' 'scan_pct' ...
        'recon_res_x' 'recon_res_y'}];
end
header.tblcols = [header.tblcols {'rescale_int' 'rescale_slope' ...
    'scale_slope' 'window_center' 'window_width' 'angulation_ap' ...
    'angulation_fh' 'angulation_rl' 'offcentre_ap' 'offcentre_fh' ...
    'offcentre_rl'}];
if (header.par_version > 4)
    header.tblcols = [header.tblcols {'thick' 'gap'}];
end
header.tblcols = [header.tblcols {'disp_orien' 'slice_orien' 'fmri' ...
    'end_dia_end_sys' 'spacing_x' 'spacing_y' 'echo_time' 'dtime' 'ttime' ...
    'diff_b_fac' }];
if (header.par_version > 4)
    header.tblcols = [header.tblcols {'number_ave'}];
end
header.tblcols = [header.tblcols {'flip' }];
if (header.par_version > 4)
    header.tblcols = [header.tblcols {'card_freq' 'min_RR_int' 'max_RR_int'...
        'turbo' 'inv_delay' 'b_value' 'grad_ori' 'contrast' 'anisotropy' ...
        'diffusion_ap' 'diffusion_fh' 'diffusion_rl'}];
end
if (header.par_version == 4.2)
    header.tblcols = [header.tblcols {'ASL_label'}];
end
header.tblcols = cell2struct (num2cell(1:length(header.tblcols)), ...
    header.tblcols, 2);

% sizes and labels for image dimensions
header.sizes = [header.nrows header.ncols header.nslices header.nechoes ...
    header.ndynamics length(unique(A(:,5))) header.nphases];
header.dims = {'rows', 'cols', 'slices', 'echoes', ...
    'dynamics' 'dataTypes' 'phases'};
if (header.par_version > 4)
    header.sizes = [header.sizes length(unique(A(:,43)))]; %grad ori
    header.dims = [header.dims {'diffusion'}];
    if (header.par_version == 4.2)
        header.sizes = [header.sizes header.nASL_type];
        header.dims = [header.dims {'aslTypes'}];
    end
end

% reorder fields in struct
header = orderfields(header);

%% display unless quiet
if ~quiet
    fprintf(['Reading in %s:\n nrows\t= %d\n ncols\t= %d\n' ...
        ' nslices\t= %d\n nechoes\t= %d\n ndynamics\t= %d\n' ...
        ' types\t= %d\n phases\t= %d\n'], ...
        filename_par, header.sizes(1:7));
    if( header.par_version >= 4)
        fprintf(' bdirs\t= %d\n', header.sizes(8));
    end
    if( header.par_version == 4.2)
        fprintf(' nASL_type\t= %d\n', header.sizes(9));
    end
end

% bottom table array elements:
% v3
% sl ec dyn ph ty(2) idx (re)scale(3) window(2) angulation(3) offcentre(3)
% info(4) spacing(2) echo dtime ttime diff flip
% v4.1
% sl ec dyn ph ty(2) idx pix scan% rec-size(2) (re)scale(3) window(2)
% angulation(3) offcentre(3) thick gap info(4) spacing(2) echo dtime ttime
% diff avg flip freq RR-int(2) turbo delay b grad cont anis diffusion(3)
% v4.2
% add L.ty (label type) to end