function varargout = fourDvis(varargin)
% FOURDVIS MATLAB code for fourDvis.fig
%      FOURDVIS, by itself, creates a new FOURDVIS or raises the existing
%      singleton*.
%
%      H = FOURDVIS returns the handle to a new FOURDVIS or the handle to
%      the existing singleton*.
%
%      FOURDVIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FOURDVIS.M with the given input arguments.
%
%      FOURDVIS('Property','Value',...) creates a new FOURDVIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fourDvis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fourDvis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fourDvis

% Last Modified by GUIDE v2.5 19-Jun-2020 11:04:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fourDvis_OpeningFcn, ...
                   'gui_OutputFcn',  @fourDvis_OutputFcn, ...
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


% --- Executes just before fourDvis is made visible.
function fourDvis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fourDvis (see VARARGIN)
global PlanesVis branchListVis segmentVis directory res caseFilePath
global pVis figVis dcm_obj_vis

PlanesVis = varargin{1}; %corners for all planes
branchListVis = varargin{2}; %locations/labels for all vessel points
segmentVis = varargin{3}; %segmentation mask
caseFilePath = varargin{4}; %directory of pcviprData file (imageData)
res = varargin{5}; %pixel resolution (mm)
pVis = []; %here to remove warning, used in myupdatefcn

% Choose default command line output for fourDvis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
movegui(handles.figure1,'southwest'); %move to bottom left (WORK)
%set(handles.figure1,'Position',[281 10 190 18]); %move to bottom left (HOME)

%%% Loading Data
load(caseFilePath,'imageData') %load data_struct

%%% Save all data in handles structure
handles.MAG = imageData.MAG;
handles.CD = imageData.CD;
handles.U = imageData.V(:,:,:,1); 
handles.V = imageData.V(:,:,:,2); 
handles.W = imageData.V(:,:,:,3); 
handles.res = size(imageData.MAG); %new matrix sizes (after cropping)
handles.MIMICS.delX = res; %pixel size for mimics masks
handles.MIMICS.delY = res;
handles.MIMICS.delZ = res;
handles.MaskID_struct(1).IDX = find(segmentVis); %index active pixels
handles.MaskID_struct(1).branchListVis = branchListVis';
handles.MaskID_struct(1).PlanesVis = PlanesVis;
handles.MaskID_struct(1).VECvis = 0; %show velocity flag
handles.MaskID_struct(1).ISOvar.vis = 1; %show isosurface flag
handles.MaskID_struct(1).ISOvar.color = 'red'; %isosurface color
handles.MaskID_struct(1).ISOvar.alpha = 1; %isosurface transparency
handles.MaskID_struct(1).ISOvar.colorval = 1; %color value (in list)
handles.MaskID_struct(1).ISOvar.smoothval = 1; %smooth value (in list)

clear imageData

%%% Create magnitude image planes for slice visualizations
figVis = figure('CloseRequestFcn',@my_closereq); %open second figure
set(figVis,'Position',[1006 43 902 692]); %WORK
%set(figVis,'Position',[2367 20 902 692]); %HOME

[X,Y,Z] = meshgrid(1:size(handles.MAG,1),1:size(handles.MAG,2),1:size(handles.MAG,3)); 
hold on     
hsurfacesx = slice(X,Y,Z,permute(handles.MAG,[2 1 3]),1,[],[]); %add planes to images
hsurfacesy = slice(X,Y,Z,permute(handles.MAG,[2 1 3]),[],1,[]);
hsurfacesz = slice(X,Y,Z,permute(handles.MAG,[2 1 3]),[],[],1);

% Turn slice visualization off at first 
set(hsurfacesx,'FaceColor','interp','EdgeColor','none','visible','off','PickableParts','none')
set(hsurfacesy,'FaceColor','interp','EdgeColor','none','visible','off','PickableParts','none')
set(hsurfacesz,'FaceColor','interp','EdgeColor','none','visible','off','PickableParts','none')
colormap gray
caxis auto
hold off
freezeColors %freeze colormap of slices

handles.hsurfacesy = hsurfacesy;
handles.hsurfacesx = hsurfacesx;
handles.hsurfacesz = hsurfacesz;

% Set slider increments to allow for vis of every slice
steps = [1/(handles.res(1)-1) 10/(handles.res(1)-1)]; %x slider
set(handles.xplaneloc,'SliderStep',steps);
set(handles.xplaneloc,'Min',0);
steps = [1/(handles.res(2)-1) 10/(handles.res(2)-1)]; %y slider
set(handles.yplaneloc,'SliderStep',steps);
set(handles.yplaneloc,'Min',0);
steps = [1/(handles.res(3)-1) 10/(handles.res(3)-1)]; %z slider
set(handles.zplaneloc,'SliderStep',steps);
set(handles.zplaneloc,'Min',0);


%%% Plot isosurface over current figure
figure(figVis)
handles.MaskID_struct(1).ISOsurf = patch(isosurface(permute(segmentVis,[2 1 3]),0.5), ... 
    'FaceAlpha',handles.MaskID_struct(1).ISOvar.alpha);
set(handles.MaskID_struct(1).ISOsurf,'FaceColor',handles.MaskID_struct(1).ISOvar.color, ...
    'EdgeColor','none','PickableParts','none');
set(gcf,'color','black'); %set background to black
axis off tight %remove axis
view([1 0.1 0.1]); %set to sag.(offset 0.1 to see cor. and ax. MAG planes)
axis vis3d
daspect([1 1 1]) %set aspect ratio so not elongated
set(gca,'zdir','reverse') %flip angiogram upside down (correct orientation)
camlight headlight; %make isosurface shine
lighting gouraud %smooth shine

[X,Y,Z] = meshgrid(1:size(handles.MAG,2),1:size(handles.MAG,1),1:size(handles.MAG,3)); 
%%% Plot velocity vectors
valsUse = logical(segmentVis); %binarize segmentation
hold on
vecLength = round(get(handles.vecLengthSlider,'Value')); %vector density
handles.q = quiver3(Y(valsUse),X(valsUse),Z(valsUse), ... 
    -handles.V(valsUse),-handles.U(valsUse),-handles.W(valsUse),vecLength);
hold off

% Compute magnitude of the vectors
handles.mags = sqrt(sum(cat(2,handles.q.UData(:),handles.q.VData(:), ...
    reshape(handles.q.WData,numel(handles.q.UData),[])).^2,2));

% Now determine the color to make each arrow using a colormap
currentColormap = colormap('jet'); %set initial colormap
[~,~,ind] = histcounts(handles.mags,size(currentColormap,1));
cmap = uint8(ind2rgb(ind(:),currentColormap)*255); %colormap to get RGB
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap,[1 3 1]),[2 1 3]);

% Repeat color 3 times (using 1:3 below) b/c each arrow has 3 vertices
set(handles.q.Head,'ColorBinding','interpolated', ... 
    'ColorData',reshape(cmap(1:3,:,:),[],4).');

% Repeat color twice (using 1:2 below) b/c each tail has 2 vertices
set(handles.q.Tail,'ColorBinding','interpolated', ...
    'ColorData',reshape(cmap(1:2,:,:),[],4).');

set(gca,'CLim',[min(handles.mags),max(handles.mags)]); %auto set min/max
set(handles.q,'PickableParts','none','visible','off') %turn off initially
handles.c = colorbar; %save colorbar in handle
handles.c.Color = [1 1 1]; %white
handles.c.LineWidth = 3; %width of colorbar
handles.c.FontSize = 20; %size of units displayed
ylabel(handles.c,'Velocity cm/sec') %colorbar caption
set(handles.c,'visible','off') %turn off initially

handles.cMIN = min(handles.mags); %save min vector magnit. as colorbar min
handles.cMAX = max(handles.mags); %save max vector magnit. as colorbar max


%%% Plot data cursor plane
hold on 
handles.hscatter = scatter3(branchListVis(:,1),branchListVis(:,2),branchListVis(:,3), ... 
    'MarkerFaceAlpha',0,'MarkerEdgeAlpha',0); 
hold off
dcm_obj_vis = datacursormode(figVis);
datacursormode on;
dcm_obj_vis.DisplayStyle = 'window';
set(dcm_obj_vis,'UpdateFcn',@myupdatefcn_planes)
dcm_obj_vis.createDatatip(handles.hscatter);
handles.StructLoc = 1;

%%% Set to old toolbar style
figVis.CurrentAxes.Toolbar = [];
addToolbarExplorationButtons(figVis)
updateDataCursors(dcm_obj_vis)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fourDvis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fourDvis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: delete(hObject) closes the figure
global figVis
delete(figVis)
delete(hObject);



% --- Executes on selection change in maskselection.
function maskselection_Callback(hObject, eventdata, handles)
% hObject    handle to maskselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns maskselection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from maskselection
handles.velonoff.Value = handles.MaskID_struct(handles.maskselection.Value).VECvis;
handles.isovisible.Value = handles.MaskID_struct(handles.maskselection.Value).ISOvar.vis;
handles.colorselection.Value = handles.MaskID_struct(handles.maskselection.Value).ISOvar.colorval;
handles.isoaplha.Value = handles.MaskID_struct(handles.maskselection.Value).ISOvar.alpha;
handles.smoothiso.Value = handles.MaskID_struct(handles.maskselection.Value).ISOvar.smoothval;
                
% --- Executes during object creation, after setting all properties.
function maskselection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addmask.
function addmask_Callback(hObject, eventdata, handles)
% hObject    handle to addmask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global figVis PlanesVis branchListVis res
figure(figVis)

%%% Read Mimics segmentation files
% User grabs text file from folder
[FileName,PathName] = uigetfile('*.txt','Select the Mimics Mask files','MultiSelect', 'on');

loopNum = 1; %default number of files
if iscell(FileName) %if more than one file is selected
    loopNum = size(FileName,2); %loopNum is now that # of files
end

for n = 1:loopNum %for each segmentation mask file
    if iscell(FileName)
        splitStr = regexp(FileName{n},'.txt','split');
        [xM,yM,zM,~] = textread([PathName FileName{n}],'%f,%f,%f,%f');
    else
        splitStr = regexp(FileName,'.txt','split');
        [xM,yM,zM,~] = textread([PathName FileName],'%f,%f,%f,%f');
    end

    badMaskFlag = max(xM)>size(handles.MAG,1) || max(yM)>size(handles.MAG,2) ...
    || max(zM)>size(handles.MAG,3);
    if badMaskFlag
        set(handles.errorText,'String','Mask is outside of image bounds. Try another mask')
        pause(3);
        set(handles.errorText,'String','');
        break
    end 
    
    Maskname = splitStr{1}; %write mask name
    contents = cellstr(handles.maskselection.String);

    if size(contents,1)<9 %if selected less than 9 files
        handles.StructLoc = size(contents,1)+1; %CD is first file
        
        % Add the name to the drop down menu
        contents{handles.StructLoc} = Maskname;
        handles.maskselection.String = contents;

        % Gets the x,y,z cordinates from a mimics mask file
        xM = xM/handles.MIMICS.delX+1;
        yM = yM/handles.MIMICS.delY+1;
        zM = zM/handles.MIMICS.delZ+1;

        idx = sub2ind(size(handles.MAG),round(xM),round(yM),round(zM));
        handles.MaskID_struct(handles.StructLoc).IDX = idx';

        %%% Create centerlines from mask
        sortingCriteria = 3;
        spurLength = 8;
        segmenttemp = zeros(size(handles.MAG));
        segmenttemp(handles.MaskID_struct(handles.StructLoc).IDX) = 1;
        vMean = cat(4,handles.U,handles.V,handles.W);
        [~,~,branchListVis,~] = feature_extraction(sortingCriteria,spurLength,vMean,segmenttemp);

        branchListSmooth = ones([size(branchListVis,1),size(branchListVis,2)]);
        for q = 1:max(branchListVis(:,4))
            branchActual = branchListVis(branchListVis(:,4)==q,:);
            if ~isempty(branchActual)
                xyz = [branchActual(:,1)';branchActual(:,2)';branchActual(:,3)'];

                % Given data in a point matrix, xyz, which is 3 x number of points
                [ndim,npts]=size(xyz);
                xyzp=zeros(size(xyz));

                % Cubic spline smoothing
                for k=1:ndim
                    pp = csaps(1:npts,xyz(k,:),0.3750);
                    xyzp(k,:)=ppval(pp,1:npts);
                end
              branchListSmooth(branchListVis(:,4)==q,1:3) = xyzp';
              branchListSmooth(branchListVis(:,4)==q,4:5) = branchListVis(branchListVis(:,4)==q,4:5);
            end
        end

        % Fills in the display metrics for all different masks
        handles.MaskID_struct(handles.StructLoc).branchListVis = branchListSmooth';
        handles.MaskID_struct(handles.StructLoc).PlanesVis = makeITPlaneGUI(branchListSmooth)';
        handles.MaskID_struct(handles.StructLoc).VECvis = 0;
        handles.MaskID_struct(handles.StructLoc).ISOvar.vis = 1;
        handles.MaskID_struct(handles.StructLoc).ISOvar.color = 'red';
        handles.MaskID_struct(handles.StructLoc).ISOvar.colorval = 1;
        handles.MaskID_struct(handles.StructLoc).ISOvar.smoothval = 1;
        handles.MaskID_struct(handles.StructLoc).ISOvar.alpha = 1;

        % Create isosurfaces
        hold on
        handles.MaskID_struct(handles.StructLoc).ISOsurf = ...
            patch(isosurface(segmenttemp,0.5),'FaceAlpha', ...
            handles.MaskID_struct(handles.StructLoc).ISOvar.alpha);
        set(handles.MaskID_struct(handles.StructLoc).ISOsurf,...
            'FaceColor',handles.MaskID_struct(handles.StructLoc).ISOvar.color, ...
            'EdgeColor', 'none','PickableParts','none');
        hold off   
    else
        disp('Can only select up to 9 masks');
    end
end

if ~badMaskFlag
    IDuse = [];
    for z = 1:handles.StructLoc
        IDuse = [IDuse;handles.MaskID_struct(z).ISOvar.vis,handles.MaskID_struct(z).VECvis];
    end
    IDuse = find(sum(IDuse,2));

    Allbranch = [handles.MaskID_struct(IDuse).branchListVis]';
    [~,IDloc] = unique(Allbranch(:,1:3),'rows');
    branchListVis = Allbranch(IDloc,:);
    Allplanes =  [handles.MaskID_struct(IDuse).PlanesVis];
    PlanesVis = Allplanes(IDloc,:);
    PlanesVis = reshape(PlanesVis,[size(PlanesVis,1),4,3]);
    set(handles.hscatter,'XData',branchListVis(:,1),'YData',branchListVis(:,2),'ZData',branchListVis(:,3));
    datacursormode on;

    delete(findall(gcf,'Type','light'))
    camlight headlight;
    lighting gouraud

    set(handles.hsurfacesx,'FaceLighting','none')
    set(handles.hsurfacesy,'FaceLighting','none')
    set(handles.hsurfacesz,'FaceLighting','none')
end 



guidata(hObject, handles);



% --- Executes on button press in isovisible.
function isovisible_Callback(hObject, eventdata, handles)
% hObject    handle to isovisible (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of isovisible
global PlanesVis branchListVis 

currMask = handles.maskselection.Value; %grab current selected mask
val = get(hObject,'Value'); %see if isosurfaces are visible
if val == 0 %if the box is unselected
    set(handles.MaskID_struct(currMask).ISOsurf ,'visible','off')
    handles.MaskID_struct(currMask).ISOvar.vis = 0;
else
    set(handles.MaskID_struct(currMask).ISOsurf,'visible','on')
    handles.MaskID_struct(currMask).ISOvar.vis = 1;
end

IDuse = [];
for z = 1:handles.StructLoc
    IDuse = [IDuse;handles.MaskID_struct(z).ISOvar.vis,handles.MaskID_struct(z).VECvis];
end
IDuse = find(sum(IDuse,2));

if ~isempty(IDuse)
    Allbranch = [handles.MaskID_struct(IDuse).branchListVis]';
    [~,IDloc] = unique(Allbranch(:,1:3),'rows');
    branchListVis = Allbranch(IDloc,:);
    Allplanes =  [handles.MaskID_struct(IDuse).PlanesVis];
    PlanesVis = Allplanes(IDloc,:);
    PlanesVis = reshape(PlanesVis,[size(PlanesVis,1),4,3]);
    set(handles.hscatter,'XData',branchListVis(:,1),'YData',branchListVis(:,2),'ZData',branchListVis(:,3));
end

guidata(hObject, handles);



% --- Executes on selection change in colorselection.
function colorselection_Callback(hObject, eventdata, handles)
% hObject    handle to colorselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns colorselection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from colorselection
contents = cellstr(get(hObject,'String'));
currentColor = contents{get(hObject,'Value')};
set(handles.MaskID_struct(handles.maskselection.Value).ISOsurf , ...
    'FaceColor',currentColor,'EdgeColor', 'none','PickableParts','none');
handles.MaskID_struct(handles.maskselection.Value).ISOvar.colorval = handles.colorselection.Value;
handles.MaskID_struct(handles.maskselection.Value).ISOvar.color = currentColor;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function colorselection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function isoaplha_Callback(hObject, eventdata, handles)
% hObject    handle to isoaplha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.MaskID_struct(handles.maskselection.Value).ISOsurf ,'FaceAlpha',get(hObject,'Value'));
handles.MaskID_struct(handles.maskselection.Value).ISOvar.alpha = handles.isoaplha.Value;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function isoaplha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isoaplha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
    set(hObject, 'Min', 0);
    set(hObject, 'Max', 1);
end



% --- Executes on selection change in smoothiso.
function smoothiso_Callback(hObject, eventdata, handles)
% hObject    handle to smoothiso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns smoothiso contents as cell array
%        contents{get(hObject,'Value')} returns selected item from smoothiso
global figVis

contents = cellstr(get(hObject,'String'));
smoothfactor = contents{get(hObject,'Value')};

figure(figVis)
delete(findall(gcf,'Type','light'))
delete(handles.MaskID_struct(handles.maskselection.Value).ISOsurf )

segmentVis = zeros(size(handles.MAG));
segmentVis(handles.MaskID_struct(handles.maskselection.Value).IDX) = 1;

if smoothfactor == '1'
    handles.MaskID_struct(handles.maskselection.Value).ISOsurf  =  ...
        patch(isosurface(permute(segmentVis,[2 1 3]),0.5),'FaceAlpha',get(handles.isoaplha, 'Value'));
    names = get(handles.colorselection, 'String');
    current = get(handles.colorselection, 'Value');
    set(handles.MaskID_struct(handles.maskselection.Value).ISOsurf , ...
        'FaceColor',names{current},'EdgeColor', 'none','PickableParts','none');
    axis off tight
    axis vis3d
    daspect([1 1 1])
    camlight headlight;
    lighting gouraud
else
    handles.MaskID_struct(handles.maskselection.Value).ISOsurf = ... 
        patch(isosurface(smooth3(permute(segmentVis,[2 1 3]),'box',str2double(smoothfactor)),0.25), ... 
        'FaceAlpha',get(handles.isoaplha, 'Value'));
    names = get(handles.colorselection, 'String');
    current = get(handles.colorselection, 'Value');
    set(handles.MaskID_struct(handles.maskselection.Value).ISOsurf , ...
        'FaceColor',names{current},'EdgeColor', 'none','PickableParts','none');
    axis off tight
    axis vis3d
    daspect([1 1 1])
    camlight headlight;
    lighting gouraud
end

handles.MaskID_struct(handles.maskselection.Value).ISOvar.smoothval = handles.smoothiso.Value;

if handles.MaskID_struct(handles.maskselection.Value).ISOvar.vis == 0 
    set(handles.MaskID_struct(handles.maskselection.Value).ISOsurf ,'visible','off')
else
    set(handles.MaskID_struct(handles.maskselection.Value).ISOsurf ,'visible','on')
end
    
set(handles.hsurfacesx,'FaceLighting','none')
set(handles.hsurfacesy,'FaceLighting','none')
set(handles.hsurfacesz,'FaceLighting','none')

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function smoothiso_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smoothiso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in velonoff.
function velonoff_Callback(hObject, eventdata, handles)
% hObject    handle to velonoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of velonoff
% Hint: get(hObject,'Value') returns toggle state of isovisible

global figVis PlanesVis branchListVis segmentVis
figure(figVis)

handles.MaskID_struct(handles.maskselection.Value).VECvis = handles.velonoff.Value;

IDuse = [];
for z = 1:handles.StructLoc
    IDuse = [IDuse;handles.MaskID_struct(z).VECvis];
end
IDuse = find(IDuse);

if isempty(IDuse)
    set(handles.q,'visible','off')
    set(handles.c,'visible','off')
else
    AllIdx = unique([handles.MaskID_struct(IDuse).IDX]);
    [Y,X,Z] = ind2sub(size(handles.MAG),AllIdx);

    %%% Color Calcs
    % Compute magnitude of the vectors
    mags = sqrt(sum(cat(2, -handles.U(AllIdx), -handles.V(AllIdx), ...
                reshape(-handles.W(AllIdx), numel(-handles.U(AllIdx)), [])).^2, 2));

    % Keeps all values
    mags(mags<handles.cMIN) = handles.cMIN;
    mags(mags>handles.cMAX) = handles.cMAX;

    contents = cellstr(get(handles.colorscalevec,'String'));
    CMAP = contents{get(handles.colorscalevec,'Value')};
    currentColormap = colormap(CMAP);

    % Now determine the color to make each arrow using a colormap
    [~, ~, ind] = histcounts(mags, size(currentColormap, 1));

    % Now map this to a colormap to get RGB
    cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
    cmap(:,:,4) = 255;
    cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

    % Repeat color 3 times (using 1:3 below) b/c each arrow has 3 vertices
    set(handles.q.Head, ...
        'ColorBinding', 'interpolated', ...
        'ColorData', reshape(cmap(1:3,:,:), [], 4).');

    % Repeat color twice (using 1:2 below) b/c each tail has 2 vertices
    set(handles.q.Tail,'ColorBinding', 'interpolated', ...
        'ColorData', reshape(cmap(1:2,:,:), [], 4).');

    %set(handles.q,'XData',X,'YData',Y,'ZData',Z,'UData',-handles.V(AllIdx), ...
        %'VData',-handles.U(AllIdx),'WData',-handles.W(AllIdx),'PickableParts','none')
    set(handles.q,'visible','on')
    set(handles.c,'visible','on')
    set(gca, 'CLim', [min(mags), max(mags)]);

    IDuse = [];
    for z = 1:handles.StructLoc
        IDuse = [IDuse;handles.MaskID_struct(z).ISOvar.vis,handles.MaskID_struct(z).VECvis];
    end
    IDuse = find(sum(IDuse,2));
    
    Allbranch = [handles.MaskID_struct(IDuse).branchListVis]';
    [~,IDloc] = unique(Allbranch(:,1:3),'rows');
    branchListVis = Allbranch(IDloc,:);
    Allplanes =  [handles.MaskID_struct(IDuse).PlanesVis];
    PlanesVis = Allplanes(IDloc,:);
    PlanesVis = reshape(PlanesVis,[size(PlanesVis,1),4,3]);
    
    hold on 
    handles.hscatter = scatter3(branchListVis(:,1),branchListVis(:,2),branchListVis(:,3), ... 
    'MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
    hold off
    
    set(handles.hscatter,'XData',branchListVis(:,1),'YData',branchListVis(:,2),'ZData',branchListVis(:,3));
end

guidata(hObject, handles);



function colorMIN_Callback(hObject, eventdata, handles)
% hObject    handle to colorMIN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of colorMIN as text
%        str2double(get(hObject,'String')) returns contents of colorMIN as a double
global figVis

contents = cellstr(get(handles.colorscalevec,'String'));
CMAP = contents{get(handles.colorscalevec,'Value')};

figure(figVis)
cMIN = str2double(get(hObject,'String'));
handles.cMIN = cMIN;
delete(handles.c)
currentColormap = colormap(CMAP);

handles.mags = sqrt(sum(cat(2, handles.q.UData(:), handles.q.VData(:), ...
            reshape(handles.q.WData, numel(handles.q.UData), [])).^2, 2));      
handles.mags(handles.mags<handles.cMIN) = handles.cMIN;
handles.mags(handles.mags>handles.cMAX) = handles.cMAX;

% Now determine the color to make each arrow using a colormap
[~, ~, ind] = histcounts(handles.mags, size(currentColormap, 1));

% Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

% Repeat color 3 times (using 1:3 below) b/c each arrow has 3 vertices
set(handles.q.Head,'ColorBinding', 'interpolated','ColorData', ...
    reshape(cmap(1:3,:,:), [], 4).');

% Repeat color twice (using 1:2 below) b/c each tail has 2 vertices
set(handles.q.Tail,'ColorBinding', 'interpolated','ColorData', ...
    reshape(cmap(1:2,:,:), [], 4).');

set(gca, 'CLim', [min(handles.mags), max(handles.mags)]);
set(handles.q,'PickableParts','none')

handles.c = colorbar;
handles.c.Color = [1 1 1];
handles.c.LineWidth = 3;
handles.c.FontSize = 20;
ylabel(handles.c, 'Velocity cm/sec')

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function colorMIN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorMIN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function colorMAX_Callback(hObject, eventdata, handles)
% hObject    handle to colorMAX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of colorMAX as text
%        str2double(get(hObject,'String')) returns contents of colorMAX as a double

global figVis

contents = cellstr(get(handles.colorscalevec,'String'));
CMAP = contents{get(handles.colorscalevec,'Value')};

figure(figVis)
cMAX = str2double(get(hObject,'String'));
handles.cMAX = cMAX;

delete(handles.c)

currentColormap = colormap(CMAP);

handles.mags = sqrt(sum(cat(2, handles.q.UData(:), handles.q.VData(:), ...
            reshape(handles.q.WData, numel(handles.q.UData), [])).^2, 2));
        
handles.mags(handles.mags<handles.cMIN) = handles.cMIN;
handles.mags(handles.mags>handles.cMAX) = handles.cMAX;

% Now determine the color to make each arrow using a colormap
[~, ~, ind] = histcounts(handles.mags, size(currentColormap, 1));

% Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

% Repeat color 3 times (using 1:3 below) b/c each arrow has 3 vertices
set(handles.q.Head,'ColorBinding', 'interpolated','ColorData', ...
    reshape(cmap(1:3,:,:), [], 4).');

% Repeat color twice (using 1:2 below) b/c each tail has 2 vertices
set(handles.q.Tail,'ColorBinding', 'interpolated','ColorData', ...
    reshape(cmap(1:2,:,:), [], 4).');

set(gca, 'CLim', [min(handles.mags), max(handles.mags)]);
set(handles.q,'PickableParts','none')

handles.c = colorbar;
handles.c.Color = [1 1 1];
handles.c.LineWidth = 3;
handles.c.FontSize = 20;
ylabel(handles.c, 'Velocity cm/sec')
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function colorMAX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorMAX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in colorscalevec.
function colorscalevec_Callback(hObject, eventdata, handles)
% hObject    handle to colorscalevec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns colorscalevec contents as cell array
%        contents{get(hObject,'Value')} returns selected item from colorscalevec
global figVis
figure(figVis)

contents = cellstr(get(hObject,'String'));
CMAP = contents{get(hObject,'Value')};
delete(handles.c)
figure(figVis)

IDuse = [];
for z = 1:handles.StructLoc
    IDuse = [IDuse;handles.MaskID_struct(z).VECvis];
end
IDuse = find(IDuse);

if ~isempty(IDuse)
    AllIdx = unique([handles.MaskID_struct(IDuse).IDX]);
    [Y,X,Z] = ind2sub(size(handles.MAG),AllIdx);

    %%% Color Calcs
    % Compute magnitude of the vectors
    mags = sqrt(sum(cat(2, -handles.U(AllIdx), -handles.V(AllIdx), ...
                reshape(-handles.W(AllIdx), numel(-handles.U(AllIdx)), [])).^2, 2));
    mags(mags<handles.cMIN) = handles.cMIN;
    mags(mags>handles.cMAX) = handles.cMAX;

    currentColormap = colormap(CMAP);

    % Now determine the color to make each arrow using a colormap
    [~,~,ind] = histcounts(mags,size(currentColormap,1));

    % Now map this to a colormap to get RGB
    cmap = uint8(ind2rgb(ind(:),currentColormap)*255);
    cmap(:,:,4) = 255;
    cmap = permute(repmat(cmap,[1 3 1]),[2 1 3]);

    % Repeat color 3 times (using 1:3 below) b/c each arrow has 3 vertices
    set(handles.q.Head, ...
        'ColorBinding','interpolated', ...
        'ColorData',reshape(cmap(1:3,:,:),[],4).');   %'

    % Repeat color twice (using 1:2 below) b/c each tail has 2 vertices
    set(handles.q.Tail,'ColorBinding','interpolated', ...
        'ColorData',reshape(cmap(1:2,:,:),[],4).');

    %set(handles.q,'XData',X,'YData',Y,'ZData',Z, ... 
        %'UData',-handles.V(AllIdx),'VData',-handles.U(AllIdx), ... 
        %'WData',-handles.W(AllIdx),'PickableParts','none')

    set(gca,'CLim',[min(handles.mags), max(handles.mags)]);
    set(handles.q,'PickableParts','none')

    handles.c = colorbar;
    handles.c.Color = [1 1 1];
    handles.c.LineWidth = 3;
    handles.c.FontSize = 20;
    ylabel(handles.c,'Velocity cm/sec')

    set(handles.q,'visible','on')
    set(handles.c,'visible','on')
end

guidata(hObject, handles);
    
    
    
% --- Executes during object creation, after setting all properties.
function colorscalevec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorscalevec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function vecLengthSlider_Callback(hObject, eventdata, handles)
% hObject    handle to vecLengthSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global figVis PlanesVis branchListVis segmentVis
figure(figVis)

%delete(handles.q);

IDuse = [];
for z = 1:handles.StructLoc
    IDuse = [IDuse;handles.MaskID_struct(z).ISOvar.vis,handles.MaskID_struct(z).VECvis];
end
IDuse = find(sum(IDuse,2));

if isempty(IDuse)
    set(handles.q,'visible','off')
    set(handles.c,'visible','off')
else
    AllIdx = unique([handles.MaskID_struct(IDuse).IDX]);
    [Y,X,Z] = ind2sub(size(handles.MAG),AllIdx);
    
    %%% Color Calcs
    % Compute magnitude of the vectors
    hold on
    vecLength = round(get(handles.vecLengthSlider,'Value'));
    handles.q = quiver3(Y,X,Z,-handles.V(AllIdx),-handles.U(AllIdx), ...
        -handles.W(AllIdx),vecLength);
    hold off

    mags = sqrt(sum(cat(2, -handles.U(AllIdx), -handles.V(AllIdx), ...
                reshape(-handles.W(AllIdx), numel(-handles.U(AllIdx)), [])).^2, 2));

    % Keeps all values
    mags(mags<handles.cMIN) = handles.cMIN;
    mags(mags>handles.cMAX) = handles.cMAX;

    contents = cellstr(get(handles.colorscalevec,'String'));
    CMAP = contents{get(handles.colorscalevec,'Value')};
    currentColormap = colormap(CMAP);

    % Now determine the color to make each arrow using a colormap
    [~, ~, ind] = histcounts(mags, size(currentColormap, 1));

    % Now map this to a colormap to get RGB
    cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
    cmap(:,:,4) = 255;
    cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

    % Repeat color 3 times (using 1:3 below) b/c each arrow has 3 vertices
    set(handles.q.Head, ...
        'ColorBinding', 'interpolated', ...
        'ColorData', reshape(cmap(1:3,:,:), [], 4).');

    % Repeat color twice (using 1:2 below) b/c each tail has 2 vertices
    set(handles.q.Tail,'ColorBinding', 'interpolated', ...
        'ColorData', reshape(cmap(1:2,:,:), [], 4).');

    %set(handles.q,'XData',X,'YData',Y,'ZData',Z,'UData',-handles.V(AllIdx), ...
        %'VData',-handles.U(AllIdx),'WData',-handles.W(AllIdx),'PickableParts','none')
    set(handles.q,'visible','on')
    set(handles.c,'visible','on')
    set(gca, 'CLim', [min(mags), max(mags)]);

    IDuse = [];
    for z = 1:handles.StructLoc
        IDuse = [IDuse;handles.MaskID_struct(z).ISOvar.vis,handles.MaskID_struct(z).VECvis];
    end
    IDuse = find(sum(IDuse,2));

    Allbranch = [handles.MaskID_struct(IDuse).branchListVis]';
    [~,IDloc] = unique(Allbranch(:,1:3),'rows');
    branchListVis = Allbranch(IDloc,:);
    Allplanes =  [handles.MaskID_struct(IDuse).PlanesVis];
    PlanesVis = Allplanes(IDloc,:);
    PlanesVis = reshape(PlanesVis,[size(PlanesVis,1),4,3]);

    hold on 
    handles.hscatter = scatter3(branchListVis(:,1),branchListVis(:,2),branchListVis(:,3), ... 
    'MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
    hold off

    set(handles.hscatter,'XData',branchListVis(:,1),'YData',branchListVis(:,2),'ZData',branchListVis(:,3));
end

if ~handles.MaskID_struct(handles.maskselection.Value).VECvis
    set(handles.q,'visible','off')
    set(handles.c,'visible','off')
end 
    
    

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function vecLengthSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vecLengthSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in xplanevis.
function xplanevis_Callback(hObject, eventdata, handles)
% hObject    handle to xplanevis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of xplanevis
val = get(hObject,'Value');
if val == 0 
    set(handles.hsurfacesx,'visible','off')
else
    set(handles.hsurfacesx,'visible','on','FaceLighting','none')
    sliceNum = 1+round( get(handles.xplaneloc,'Value').*(handles.res(1)-1) );
    tempIM = squeeze(handles.MAG(sliceNum,:,:));
    NORMgary = tempIM./max(tempIM(:));
    handles.hsurfacesx.CData = cat(3, NORMgary, NORMgary, NORMgary);
    %ceil(handles.xplaneloc.Value*handles.res(1)
    handles.hsurfacesx.XData = ones(size(handles.hsurfacesx.XData)).*sliceNum;
end

guidata(hObject, handles);


 
% --- Executes on slider movement.
function xplaneloc_Callback(hObject, eventdata, handles)
% hObject    handle to xplaneloc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliceNum = 1+round( get(handles.xplaneloc,'Value').*(handles.res(1)-1) );
tempIM = squeeze(handles.MAG(sliceNum,:,:));
NORMgary = tempIM./max(tempIM(:));
handles.hsurfacesx.CData = cat(3, NORMgary, NORMgary, NORMgary);
handles.hsurfacesx.XData = ones(size(handles.hsurfacesx.XData)).*sliceNum;

% --- Executes during object creation, after setting all properties.
function xplaneloc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xplaneloc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in yplanevis.
function yplanevis_Callback(hObject, eventdata, handles)
% hObject    handle to yplanevis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of yplanevis
val = get(hObject,'Value');
if val == 0 
    set(handles.hsurfacesy,'visible','off')
else
    set(handles.hsurfacesy,'visible','on','FaceLighting','none')
    sliceNum = 1+round( get(handles.yplaneloc,'Value').*(handles.res(2)-1) );
    tempIM = squeeze(handles.MAG(:,sliceNum,:));
    NORMgary = tempIM./max(tempIM(:));
    handles.hsurfacesy.CData = cat(3, NORMgary, NORMgary, NORMgary);
    handles.hsurfacesy.YData = ones(size(handles.hsurfacesy.YData)).*sliceNum;
end

guidata(hObject, handles);


 
% --- Executes on slider movement.
function yplaneloc_Callback(hObject, eventdata, handles)
% hObject    handle to yplaneloc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliceNum = 1+round( get(handles.yplaneloc,'Value').*(handles.res(2)-1) );
tempIM = squeeze(handles.MAG(:,sliceNum,:));
NORMgary = tempIM./max(tempIM(:));
handles.hsurfacesy.CData = cat(3, NORMgary, NORMgary, NORMgary);
handles.hsurfacesy.YData = ones(size(handles.hsurfacesy.YData)).*sliceNum;

% --- Executes during object creation, after setting all properties.
function yplaneloc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yplaneloc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in zplanevis.
function zplanevis_Callback(hObject, eventdata, handles)
% hObject    handle to zplanevis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of zplanevis
val = get(hObject,'Value');
if val == 0 
    set(handles.hsurfacesz,'visible','off')
else
    set(handles.hsurfacesz,'visible','on','FaceLighting','none')
    sliceNum = 1+round( get(handles.zplaneloc,'Value').*(handles.res(3)-1) );
    tempIM = squeeze(handles.MAG(:,:,sliceNum))';
    NORMgary = tempIM./max(tempIM(:));
    handles.hsurfacesz.CData = cat(3, NORMgary, NORMgary, NORMgary);
    handles.hsurfacesz.ZData = ones(size(handles.hsurfacesz.ZData)).*sliceNum;
end

guidata(hObject, handles);



% --- Executes on slider movement.
function zplaneloc_Callback(hObject, eventdata, handles)
% hObject    handle to zplaneloc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider    
sliceNum = 1+round( get(handles.zplaneloc,'Value').*(handles.res(3)-1) );
tempIM = squeeze(handles.MAG(:,:,sliceNum))';
NORMgary = tempIM./max(tempIM(:)); %https://normangary.com/
handles.hsurfacesz.CData = cat(3, NORMgary, NORMgary, NORMgary);
handles.hsurfacesz.ZData = ones(size(handles.hsurfacesz.ZData)).*sliceNum;

% --- Executes during object creation, after setting all properties.
function zplaneloc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zplaneloc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in savefigpush.
function savefigpush_Callback(hObject, eventdata, handles)
% hObject    handle to savefigpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global figVis directory
figVis.Color = 'black';
figVis.InvertHardcopy = 'off';

if isempty(handles.savefigname.String)
    saveas(figVis,[directory '\Figure1.png']);
    set(handles.errorText,'String','Saved as Figure1.png');
    pause(3);
    set(handles.errorText,'String','');
else
   saveas(figVis,[directory '\' handles.savefigname.String '.png']);
   set(handles.errorText,'String',['Saved as ' handles.savefigname.String '.png']);
   pause(3);
   set(handles.errorText,'String','');
end



function savefigname_Callback(hObject, eventdata, handles)
% hObject    handle to savefigname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of savefigname as text
%        str2double(get(hObject,'String')) returns contents of savefigname as a double

% --- Executes during object creation, after setting all properties.
function savefigname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savefigname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%% My Update Function
function txt = myupdatefcn_planes(empt,event_obj)
% Customizes text of data tips
global PlanesVis branchListVis pVis dcm_obj_vis figVis

delete(pVis)
info_struct = getCursorInfo(dcm_obj_vis);
ptList = [info_struct.Position];
ptList = reshape(ptList,[3,numel(ptList)/3])';
pindex = zeros(size(ptList,1),1);

% Find cursor point in branchList
for n = 1:size(ptList,1)
    xIdx = find(branchListVis(:,1) == ptList(n,1));
    yIdx = find(branchListVis(xIdx,2) == ptList(n,2));
    zIdx = find(branchListVis(xIdx(yIdx),3) == ptList(n,3));
    pindex(n) = xIdx(yIdx(zIdx));
end

% Get associated branch number of full branch
bnum = branchListVis(pindex,4);

hold on  
%Update Planes for points
pVis = fill3(PlanesVis(pindex,:,1)',PlanesVis(pindex,:,2)',PlanesVis(pindex,:,3)', ...
    [0 0 1],'EdgeColor',[0 0 1],'FaceAlpha',1,'Visible','on', ...
    'PickableParts','none','Parent', figVis.Children(2)); 
% fill3(pty',ptx',ptz','r') when used with isosurface
hold off


% Get branch indices and current label point
branchActual = branchListVis(branchListVis(:,4) == bnum,5);
CurrentNum = find(branchListVis(pindex,5)==branchActual)-1;

% Update cursor text
txt = {['Branch Point: ',sprintf('%i',CurrentNum),newline ...
    'Current Branch: ',sprintf('%i',bnum)]};



%%% Close request
function my_closereq(src,callbackdata)
% Close request function 
% to display a question dialog box 
global figVis

selection = questdlg('Close This Figure?','Close Request Function',...
  'Yes','No','Yes'); 
switch selection 
  case 'Yes'
     delete(figVis)
  case 'No'
  return 
end



%%% Make Planes for centerline points
function [ Planes ] = makeITPlaneGUI( branchList )
%makeITPlane Creates a nice plane associated with the cursor point
%   Show the tangent planes that are used in the calculations for the
%   flow,area,PI,etc.

% Created by: Carson Hoffman
% Date: 03/21/2017
% University of Wisconsin Madison


%%% Getting the tangent vectors at all points
d = 2; %dist. behind/ahead of current pt for tangent plane calc (d=2->5pts)
Tangent_V = zeros(0,3); % Initialize tangent vector list
for n = 1:max(branchList(:,4))
    branchActual = branchList(branchList(:,4)==n,:);
    dir_temp = zeros(size(branchActual,1),3);
    for i = 1:size(branchActual,1)
        % Extract normal to cross-section
        if i < d+1 %if near 1st endpoint
            dir = (branchActual(i+d,1:3) - branchActual(i,1:3));
        elseif i >= size(branchActual,1)-d %if near 2nd endpoint 
            dir = (branchActual(i,1:3) - branchActual(i-d,1:3));
        else %calculate tangent from d points ahead/behind curr point
            dir = (branchActual(i+d,1:3) - branchActual(i-d,1:3));
        end
        dir_temp(i,:) = dir/norm(dir); %tangent vector with magnitude of 1
    end
    Tangent_V = [Tangent_V;dir_temp]; %add all tangents to large list
end

%%% Find normalized vector perpendicular to tangent vector
[~,idx_max] = max(abs(Tangent_V),[],2); %get max unit along rows
idx_max(idx_max==2) = 1; %flatten to 2D
max_pts = sub2ind(size(Tangent_V),[1:size(Tangent_V,1)]',idx_max);
temp = zeros(size(Tangent_V)); 
temp(max_pts) = 1; %binary matrix of location of max unit vectors
[~,idx_shift] = max(abs(circshift(temp,1,2)),[],2); %rotate (ie x->y,z->x)
shift_pts = sub2ind(size(Tangent_V),[1:size(Tangent_V,1)]',idx_shift);
V2 = zeros(size(Tangent_V));
V2(max_pts) = Tangent_V(shift_pts);
V2(shift_pts) = -Tangent_V(max_pts); % Vector 1 that is used to created the perdendicular plane
N = repmat(sqrt(sum(abs(V2).^2,2)),[1 3]); %repeat vel. magnitude as Nx3
V2 = V2./N;
V3 = cross(Tangent_V,V2);% Vector 2 that is used created the perdendicular plane
% V3,V2,Tangent_V are all orthogonal (i.e. dot( V3(1,:),Tangent_V(1,:) )=0)


%%% Interpolate
% Get the full tangent plane for all the points
r = 10; %size of plane to select from non interpolated data is r*2+1
InterpVals = 4; %choose the interpolation between points
Side = r*InterpVals; %creates correct number of points for interpolation
width = Side.*2+1; %width of plane in pixels
Mid = zeros(length(branchList),1);

% Find x values on line
temp = repmat(V2(:,1)./InterpVals,[1 Side]);
temp = cumsum(temp,2); %runs from 0 to +(r*interpVals) by unit dist/interp
temp2 = -fliplr(temp); %runs from -(r*interpVals) to 0 by unit dist/interp
x_val = [temp2 Mid temp]; %combine temps--size = N x (r*interpVals*2)+1
x_val = bsxfun(@plus,x_val,branchList(:,1)); %pointwise addition
x_val = reshape(x_val,[numel(x_val) 1]); %stretch into vector

% Find y values on line
temp = repmat(V2(:,2)./InterpVals,[1 Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
y_val = [temp2 Mid temp];
y_val = bsxfun(@plus,y_val,branchList(:,2));
y_val = reshape(y_val,[numel(y_val) 1]);

% Find z values on the line
temp = repmat(V2(:,3)./InterpVals,[1 Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
z_val = [temp2 Mid temp];
z_val = bsxfun(@plus,z_val,branchList(:,3));
z_val = reshape(z_val,[numel(z_val) 1]);

% At this point x,y,z values have created a tangent line perpendicular to
% the normal vector for all centerline points.
% Now, we begin filling out the other perpendicular line to create a plane.

% Find x values on plane
Mid = zeros(length(branchList)*(width),1);
temp = repmat(V3(:,1)./InterpVals,[width Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
x_full = [temp2 Mid temp];
x_full = bsxfun(@plus,x_full,x_val);
x_full = reshape(x_full,[length(branchList)*(width).^2,1]);

% Find y values on plane
temp = repmat(V3(:,2)./InterpVals,[(width) Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
y_full = [temp2 Mid temp];
y_full = bsxfun(@plus,y_full,y_val);
y_full = reshape(y_full,[length(branchList)*(width).^2,1]);

% Find z values on plane
temp = repmat(V3(:,3)./InterpVals,[(width) Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
z_full = [temp2 Mid temp];
z_full = bsxfun(@plus,z_full,z_val);
z_full = reshape(z_full,[length(branchList)*(width).^2,1]);

% Typecast to single and reshape
x_full = reshape(single(x_full),[length(branchList),(width).^2]);
y_full = reshape(single(y_full),[length(branchList),(width).^2]);
z_full = reshape(single(z_full),[length(branchList),(width).^2]);

% Get corners of UNINTERPOLATED planes
Planes = zeros(size(branchList,1),4,3);
Planes(:,:,1) = [x_full(:,1),x_full(:,width-InterpVals),x_full(:,end),x_full(:,end-width+1)];
Planes(:,:,2) = [y_full(:,1),y_full(:,width-InterpVals),y_full(:,end),y_full(:,end-width+1)];
Planes(:,:,3) = [z_full(:,1),z_full(:,width-InterpVals),z_full(:,end),z_full(:,end-width+1)];

