function varargout = pwvTool(varargin)
% PWVTOOL MATLAB code for pwvTool.fig
%      PWVTOOL, by itself, creates a new PWVTOOL or raises the existing
%      singleton*.
%
%      H = PWVTOOL returns the handle to a new PWVTOOL or the handle to
%      the existing singleton*.
%
%      PWVTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PWVTOOL.M with the given input arguments.
%
%      PWVTOOL('Property','Value',...) creates a new PWVTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pwvTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pwvTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pwvTool

% Last Modified by GUIDE v2.5 03-Nov-2021 20:49:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pwvTool_OpeningFcn, ...
                   'gui_OutputFcn',  @pwvTool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before pwvTool is made visible.
function pwvTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pwvTool (see VARARGIN)
global branchList res timeres flowPulsatile_val area_val 

branchList = varargin{1}; %corners for all planes
res = varargin{2}; %locations/labels for all vessel points
timeres = varargin{3}; %segmentation mask
flowPulsatile_val = varargin{4}; %directory of pcviprData file (imageData)
area_val = varargin{5}; %pixel resolution (mm)


% Choose default command line output for pwvTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pwvTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pwvTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function branchNumber_Callback(hObject, eventdata, handles)
% hObject    handle to branchNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of branchNumber as text
%        str2double(get(hObject,'String')) returns contents of branchNumber as a double


% --- Executes during object creation, after setting all properties.
function branchNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to branchNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function firstCLpoint_Callback(hObject, eventdata, handles)
% hObject    handle to firstCLpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of firstCLpoint as text
%        str2double(get(hObject,'String')) returns contents of firstCLpoint as a double


% --- Executes during object creation, after setting all properties.
function firstCLpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to firstCLpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function lastCLpoint_Callback(hObject, eventdata, handles)
% hObject    handle to lastCLpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lastCLpoint as text
%        str2double(get(hObject,'String')) returns contents of lastCLpoint as a double


% --- Executes during object creation, after setting all properties.
function lastCLpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lastCLpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in analyzeEndpointsButton.
function analyzeEndpointsButton_Callback(hObject, eventdata, handles)
% hObject    handle to analyzeEndpointsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global branchList timeres flowPulsatile_val area_val 

%text = get(handles.branchNumber,'String');
%branchNumbers = sscanf(text, '%g,');
%firstCLpoint = str2double(get(handles.firstCLpoint,'String'))+1;
%lastCLpoint = str2double(get(handles.lastCLpoint,'String'))+1;

text = get(handles.branchNumber,'String');
branchNumbers = sscanf(text, '%g,');
firstCLpoint = sscanf(get(handles.firstCLpoint,'String'),'%g,');
lastCLpoint = sscanf(get(handles.lastCLpoint,'String'),'%g,');

% Define orientation
%start = find((branchList(:,4)==branchNumbers(1)&(branchList(:,5)==firstCLpoint)));
%stop = find((branchList(:,4)==branchNumbers(end)&(branchList(:,5)==lastCLpoint)));
%if branchList(start,3)>branchList(stop,3)
%    orientSI = 1;
%else
%    orientSI = 0;
%end 

% Get flow along segments
vessel = [];
flow = [];
area = [];
for b=1:length(branchNumbers)

    % Define orientation
    start = find((branchList(:,4)==branchNumbers(b)&(branchList(:,5)==firstCLpoint(b))));
    stop = find((branchList(:,4)==branchNumbers(b)&(branchList(:,5)==lastCLpoint(b))));

    if branchList(start,3)>branchList(stop,3)
        orientSI = 1;
    else
        orientSI = 0;
    end 

    vesselInds = find(branchList(:,4)==branchNumbers(b));
    vessel_seg = branchList(vesselInds,:);
    flow_seg = flowPulsatile_val(vesselInds,:);
    area_seg = area_val(vesselInds,:);

    if orientSI
        if vessel_seg(1,3)<vessel_seg(end,3)
            vessel_seg = flipud(vessel_seg);
            flow_seg = flipud(flow_seg);
            area_seg = flipud(area_seg);
        end 
    else
        if vessel_seg(1,3)>vessel_seg(end,3)
            vessel_seg = flipud(vessel_seg);
            flow_seg = flipud(flow_seg);
            area_seg = flipud(area_seg);
        end 
    end 

    firstPointIdx = find(vessel_seg(:,5)==firstCLpoint(b));
    lastPointIdx = find(vessel_seg(:,5)==lastCLpoint(b));

    vessel = [vessel; vessel_seg(firstPointIdx:lastPointIdx,:)];
    flow =   [flow; flow_seg(firstPointIdx:lastPointIdx,:)];
    area =   [area; area_seg(firstPointIdx:lastPointIdx,:)];

end 

% Distance
positions = vessel(:,1:3);
distances = cumsum(vecnorm(diff(positions),2,2));

% maximun likelihood estimator
nFrames = size(flow,2);
scale = 5;
smoothLevel = 1; % Needs optimization ? Lets avoid further smoothing for the moment (set to 1), but 15 reduces result variability (forced?)
nFrames_interp = scale*nFrames;
timeresInt = timeres/scale;

% interp params
x = 1:nFrames;
xq = linspace(1,nFrames,nFrames_interp);

waveforms = [];
scaling = [];

for i = 1:length(flow(:,1))
    wave_smooth = smoothdata(flow(i,:),'gaussian',smoothLevel);
    %remove temporal mean, normalized to std,( to have zero mean and unit std) and interpolate
    wave_interp = interp1(x,wave_smooth,xq,'spline');
    wave_demean = wave_interp - mean(wave_interp);
    wave_norm   = wave_demean ./ std(wave_demean);
    waveforms(i,:) = wave_norm;
    
    % define scaling factor 
    wave_std = std(wave_demean);
    wave_area = area(i);
    scaling(i) = wave_area./(wave_std.^2);

end

% scaling can have issues, this just normalized outliers
scaling_th = mean(scaling) + 4*std(scaling);
scaling(scaling > scaling_th) = mean(scaling);

% scale = ones(size(distances))';
pwv_init = 5;
inParams= [waveforms(1,:)*scaling(1), pwv_init]; %[m/s]
fun1=@(inParams)PWVest3_share(inParams,distances*1e-3,waveforms(2:end,:),timeresInt*1e-3,scaling(2:end)');       
options = optimset('Display','iter', 'TolCon', 1e-7, 'TolX', 1e-7, 'TolFun', 1e-7,'DiffMinChange', 1e-3);
[params,exitflag,output] = fminunc(fun1,inParams, options);

PWV = params(end);
disp(['Maximun Likelihood Estimator from Flow: ' num2str(PWV) ' m/s']);


% --- Executes on button press in analyzeFullButton.
function analyzeFullButton_Callback(hObject, eventdata, handles)
% hObject    handle to analyzeFullButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox('Use the analyze endpoints button instead')

