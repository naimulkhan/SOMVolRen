


function varargout = SOFMFig(varargin)
% SOFMFIG MATLAB code for SOFMFig.fig
%      SOFMFIG, by itself, creates a new SOFMFIG or raises the existing
%      singleton*.
%
%      H = SOFMFIG returns the handle to a new SOFMFIG or the handle to
%      the existing singleton*.
%
%      SOFMFIG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SOFMFIG.M with the given input arguments.
%
%      SOFMFIG('Property','Value',...) creates a new SOFMFIG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SOFMFig_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SOFMFig_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SOFMFig

% Last Modified by GUIDE v2.5 26-Mar-2013 00:37:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SOFMFig_OpeningFcn, ...
    'gui_OutputFcn',  @SOFMFig_OutputFcn, ...
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


% --- Executes just before SOFMFig is made visible.
function SOFMFig_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SOFMFig (see VARARGIN)

% Choose default command line output for SOFMFig

handles.output = hObject;
handles.curdata='Engine';
handles.curdir='../MyVTKData/';
fid = fopen(strcat(handles.curdir,'filename.txt'), 'w','l');
fprintf(fid,'%s\n',handles.curdata);
fclose(fid);
load(strcat(handles.curdir,handles.curdata,'sofm.mat'));
load c3-12.mat;
load(strcat(handles.curdir,handles.curdata,'mask.mat'));
load(strcat(handles.curdir,handles.curdata,'features.mat'));

if(strcmp(handles.curdata,'INCISIX') || strcmp(handles.curdata,'INCISIX') || strcmp(handles.curdata,'CEREBRIX'))
    V=double(dicom_read_volume(strcat(handles.curdir,handles.curdata,'/','IM-0001-0001.dcm')));

else
    V=double(mha_read_volume(strcat(handles.curdir,handles.curdata,'.mha')));
end
handles.X=X;
handles.C=C;
handles.w=w;
handles.r=r;
handles.assoc=assoc;
handles.freq=freq;
handles.J=J; %the background mask, here false means foreground, see the regionGrowing folder for background sub code
handles.V=V;
handles.P=P;
handles.alreadySelected=zeros(size(X,1),1);
handles.alreadySelected=uint16(handles.alreadySelected);
handles.totClusters=1;
handles.hrange=[0.0 0.5];
handles.alphaoffset=zeros(1,5);
%set(handles.hoffsetSlider, 'Max', 1.0-handles.hrange(2));
%set(handles.hoffsetSlider, 'SliderStep', [0.05 0.05]/(1.0-handles.hrange(2)));
set(handles.checkbox1,'Value',1);
set(handles.checkbox2,'Value',1);
set(handles.checkbox3,'Value',1);
set(handles.checkbox4,'Value',1);
set(handles.checkbox5,'Value',1);

handles.selectCount=false;
handles.isclusteron=true(1,5);
handles.selectallPressed=false;
% Update handles structure
guidata(hObject, handles);
%glyphHard(X);
%axes(handles.glyphfig);
% [r,c]=stdrc(w);
% glyph(X,ones(size(X,1),1),r);
glyph(X,ones(size(X,1),1),calcColorParams(w,C,1));
%glyphHard(X,ones(size(X,1),1),w(:,2));

writeRGBAfordeselection(~handles.alreadySelected,handles);


%set(findobj('Tag','sliceSlider'),'Value',1,'Max',size(handles.V,3),'Min',1);
%set(findobj('Tag','sliceSlider'),'SliderStep',[1/(size(handles.V,3)-1) 0.1]);
% UIWAIT makes SOFMFig wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SOFMFig_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selectRegion.
function selectRegion_Callback(hObject, eventdata, handles)
% hObject    handle to selectRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




%axes(handles.glyphfig);
ss=selectGlyphPoints(handles.X,handles.alreadySelected,handles.totClusters);


handles.alreadySelected=handles.alreadySelected+ss;
if handles.selectCount
    unplot;
end
plotSelectedPoints(handles.X,handles.alreadySelected,handles.totClusters);
handles.selectCount=true;
guidata(hObject, handles);

writeRGBAforselection(ss,handles);
%DrawSlice(hObject,handles);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over selectRegion.
function selectRegion_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to selectRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function glyphfig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to glyphfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate glyphfig


% --- Executes on button press in deselectRegion.
function deselectRegion_Callback(hObject, eventdata, handles)
% hObject    handle to deselectRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%axes(handles.glyphfig);
ss=deselectGlyphPoints(handles.X,handles.alreadySelected,handles.totClusters);

handles.alreadySelected=handles.alreadySelected-ss;

if handles.selectCount
    unplot;
end

plotSelectedPoints(handles.X,handles.alreadySelected,handles.totClusters);
handles.selectCount=true;
guidata(hObject, handles);
writeRGBAfordeselection(ss,handles);






% --- Executes on button press in clearbutton.
function clearbutton_Callback(hObject, eventdata, handles)
% hObject    handle to clearbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%axes(handles.glyphfig);

handles.alreadySelected=uint16(zeros(size(handles.X,1),1));

if handles.selectCount
    for i=1:handles.totClusters
        
        curbox = findobj('Tag',strcat('checkbox',num2str(i)));
        curbutton1 = findobj('Tag',strcat('cplus',num2str(i)));
        curbutton2 = findobj('Tag',strcat('cminus',num2str(i)));
        
        set(curbox,'visible','off');
        set(curbutton1,'visible','off');
        set(curbutton2,'visible','off');
        unplot;
    end
end
handles.totClusters=1;
guidata(hObject, handles);
writeRGBAfordeselection(~handles.alreadySelected,handles);

% --- Executes on button press in selectall.
function selectall_Callback(hObject, eventdata, handles)
% hObject    handle to selectall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.selectallPressed=true;
handles.alreadySelected(~handles.alreadySelected)=uint16(handles.totClusters*ones(length(handles.alreadySelected(~handles.alreadySelected)),1));
hold on;
if handles.selectCount
    unplot;
end
hold off;

plotSelectedPoints(handles.X,handles.alreadySelected,handles.totClusters);

handles.selectCount=true;
guidata(hObject, handles);
writeRGBAforselection(handles.alreadySelected,handles);

% --- Executes on button press in genRGBA.
function genRGBA_Callback(hObject, eventdata, handles)
% hObject    handle to genRGBA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Current approach: use distance from winning node as opacity for
%selected region, use very low value(20?) for unselected region, use 0 opacity
%for background voxels


V=handles.V(:);
J=handles.J(:);
RGBA=zeros(size(V,1),4);
% j=1;

RGBA(J,:)=0;
RGBA(~J,1:3)=[handles.P(:,1) handles.P(:,2) handles.P(:,2)];
RGBA(~J,4)=50;

% for i=1:size(V,1)
%
%     if J(i,1)
%         RGBA(i,1:3)=zeros(1,3);
%         RGBA(i,4)=0;
%     end
%     if ~J(i,1)
%         RGBA(i,4)=50;
%         RGBA(i,1:3)=handles.P(j,:);
%         %RGBA(i,1:3)=handles.w(handles.assoc(j),:);
%
%         %RGBA(i,1:3)=handles.X(handles.assoc(j),:);
%         j=j+1;
%     end
%
% end

% blabla=RGBA(:,2);
% save('test.mat','blabla');
%Engine
% RGBA(:,1:3)=[(49-(RGBA(:,2)/max(RGBA(:,2)))*48) ((RGBA(:,1)/max(RGBA(:,1)))*200+2) ((RGBA(:,2)/max(RGBA(:,2)))*100+105)];

%Foot
%RGBA(:,1:3)=[(250-(RGBA(:,2)/max(RGBA(:,2)))*48) ((RGBA(:,1)/max(RGBA(:,1)))*73+178) ((RGBA(:,2)/max(RGBA(:,2)))*65+64)];

% RGBA(:,1:3)=[((RGBA(:,2)/max(RGBA(:,2)))*10+170) ((RGBA(:,1)/max(RGBA(:,1)))*10+170) ((RGBA(:,2)/max(RGBA(:,2)))*10+170)];

%Pig
RGBA(:,1:3)=[(240-(RGBA(:,2)/max(RGBA(:,2)))*60) ((RGBA(:,1)/max(RGBA(:,1)))*58+53) ((RGBA(:,2)/max(RGBA(:,2)))*80+127)];

%RGBA(:,1:3)=[RGBA(:,1)/255 RGBA(:,2)/255 RGBA(:,3)/255];
%RGBA(:,1:3)=hsv2rgb(RGBA(:,1:3)).*255;
% RGBA(:,1:3)=[(RGBA(:,1)/max(RGBA(:,1)))*61+34 (RGBA(:,2)/max(RGBA(:,2)))*107+111 (RGBA(:,3)/max(RGBA(:,3)))*64+37];

%RGBA(:,1:3)=[(RGBA(:,1)/max(RGBA(:,1)))*101+140 (RGBA(:,2)/max(RGBA(:,2)))*82+111 (RGBA(:,3)/max(RGBA(:,3)))*22+9];
%RGBA(RGBA(:,2)>0.8,1)=0;
%RGBA(:,1:3)=[((RGBA(:,1)-min(RGBA(:,1)))/(max(RGBA(:,1))-min(RGBA(:,1))))*210 ((RGBA(:,2)-min(RGBA(:,2)))/(max(RGBA(:,2))-min(RGBA(:,2))))*230 ((RGBA(:,3)-min(RGBA(:,3)))/(max(RGBA(:,3))-min(RGBA(:,1))))*200];




fid = fopen(strcat(handles.curdir,handles.curdata,'rgba.bin'), 'w','l');
fwrite(fid, RGBA', 'float','l');
fclose(fid);
uiwait(msgbox('Done!'));
%clear all;



function writeRGBAforselection(ss,handles)
assocmask=ismember(handles.assoc,find(ss)); %assoc contains the association of each voxel to a node (1..642), find(ss) returns the index of who were selected, so this line returns logical 1 for voxels selected

%sum(assocmask)
selectInd=find(~handles.J(:));
selectInd=selectInd(assocmask);
siz=size(handles.V);
[Ix,Iy,Iz]=ind2sub(siz,selectInd);
alphaval=handles.w(handles.assoc(assocmask),:)-handles.P(assocmask,:);

alphaval=sqrt(sum(abs(alphaval).^2,2));
alphaval=(1-alphaval/max(alphaval))*155+100;

% RGBA(:,1:3)=handles.P(assocmask,:);
% RGBA(:,1:3)=[(RGBA(:,2)/max(RGBA(:,2)))*11+19 (RGBA(:,1)/max(RGBA(:,1)))*81+134 (RGBA(:,3)/max(RGBA(:,3)))*60+99];
% RGBA=horzcat(Ix,Iy,Iz,RGBA,alphaval);
%
writeMat=horzcat(Ix,Iy,Iz,alphaval);

fid = fopen('../MyVTKData/select.bin', 'w','l');
fwrite(fid,0,'uint','l'); %select or cluster?
fwrite(fid, writeMat', 'float','l');
fclose(fid);
dummyval=1;
save('../MyVTKData/check.txt','dummyval','-ASCII');
while ~exist('../MyVTKData/doneren.txt','file')
    %pause(0.5);
end
delete('../MyVTKData/doneren.txt');
if exist('../MyVTKData/check.txt','file')
    delete('../MyVTKData/check.txt');
end
%delete('../MyVTKData/select.bin');

function writeRGBAfordeselection(ss,handles)
assocmask=ismember(handles.assoc,find(ss));
%sum(assocmask)
selectInd=find(~handles.J(:));
selectInd=selectInd(assocmask);
siz=size(handles.V);
[Ix,Iy,Iz]=ind2sub(siz,selectInd);
alphaval=1*ones(size(Ix,1),1);
% RGBA(:,1:3)=handles.P(assocmask,:);
% RGBA(:,1:3)=[(RGBA(:,1)/max(RGBA(:,1)))*101+140 (RGBA(:,2)/max(RGBA(:,2)))*82+111 (RGBA(:,3)/max(RGBA(:,3)))*22+9];
% RGBA=horzcat(Ix,Iy,Iz,RGBA,alphaval);
writeMat=horzcat(Ix,Iy,Iz,alphaval);

fid = fopen('../MyVTKData/select.bin', 'w','l');
fwrite(fid,0,'uint','l'); %select or cluster?
fwrite(fid,writeMat', 'float','l');
fclose(fid);
dummyval=1;
save('../MyVTKData/check.txt','dummyval','-ASCII');
while ~exist('../MyVTKData/doneren.txt','file')
    %pause(0.5);
end
delete('../MyVTKData/doneren.txt');
if exist('../MyVTKData/check.txt','file')
    delete('../MyVTKData/check.txt');
end

%delete('../MyVTKData/select.bin');
% --- Executes on button press in gencluster.
function gencluster_Callback(hObject, eventdata, handles)
% hObject    handle to gencluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for i=1:handles.totClusters
    
    curbox = findobj('Tag',strcat('checkbox',num2str(i)));
    curbutton1 = findobj('Tag',strcat('cplus',num2str(i)));
    curbutton2 = findobj('Tag',strcat('cminus',num2str(i)));
    
    set(curbox,'visible','on');
    set(curbutton1,'visible','on');
    set(curbutton2,'visible','on');
    
end
handles.totClusters=handles.totClusters+1;
guidata(hObject, handles);
calcparamsforClusters(hObject, eventdata,handles);



function calcparamsforClusters(hObject, eventdata, handles)

siz=size(handles.V);
centroid=[mean(1:siz(1)) mean(1:siz(2)) mean(1:siz(3))];
for i=1:handles.totClusters-1 %increment before calling this function!
    
    assocmask=ismember(handles.assoc,find(handles.alreadySelected==i)); %assoc contains the association of each voxel to a node (1..642), find(ss) returns the index of who were selected, so this line returns logical 1 for voxels selected
    
    %sum(assocmask)
    selectInd=find(~handles.J(:));
    selectInd=selectInd(assocmask);
    
    [Ix,Iy,Iz]=ind2sub(siz,selectInd);
    varCluster(i)=var(Ix)+var(Iy)+var(Iz);
    curcentroid=[mean(Ix) mean(Iy) mean(Iz)];
    centdist=sqrt(sum(abs(centroid-curcentroid).^2));
    handles.sat(i)=1/(1+varCluster(i));
    handles.lum(i)=1/(1+centdist);
    
end
% handles.sat=(handles.sat-min(handles.sat))./(max(handles.sat-min(handles.sat)));
% handles.lum=(handles.lum-min(handles.lum))./(max(handles.lum-min(handles.lum)));

handles.sat=handles.sat/max(handles.sat);
handles.lum=handles.lum/max(handles.lum);
for i=1:handles.totClusters-1
    if handles.sat(i)==1
        handles.sat(i)=0.8;
    end
    if handles.lum(i)==1
        handles.lum(i)=0.9;
    end
    
    handles.origopac(i)=1.2-(varCluster(i)/max(varCluster));
end
%handles.origopac
handles.huevals=linspace(handles.hrange(1),handles.hrange(2),handles.totClusters-1);
[~,handles.huevals]=deal(fix(handles.huevals), handles.huevals-fix(handles.huevals));
handles.huevals
guidata(hObject, handles);
writeHSVforClusters(handles);


function writeHSVforClusters(handles)

for i=1:handles.totClusters-1
    
    if handles.isclusteron(i)
        writeHSVforSingleCluster(handles,i);
    end
end



function writeHSVforSingleCluster(handles,i)
siz=size(handles.V);
assocmask=ismember(handles.assoc,find(handles.alreadySelected==i)); %assoc contains the association of each voxel to a node (1..642), find(ss) returns the index of who were selected, so this line returns logical 1 for voxels selected

%sum(assocmask)
selectInd=find(~handles.J(:));
selectInd=selectInd(assocmask);
%handles.huevals(i)
%handles.sat(i)
%handles.lum(i)
[Ix,Iy,Iz]=ind2sub(siz,selectInd);
RGBA=zeros(length(selectInd),4);
RGBA(:,1:3)=repmat(hsv2rgb([handles.huevals(i) handles.sat(i) handles.lum(i)]),length(selectInd),1);

% offsetvals=RGBA(1,1:3)./3;
offsetvals=RGBA(1,1:3)./4;

%RGBA(:,1:3)=[offsetvals(1)*(handles.P(assocmask,2)/max(handles.P(assocmask,2)))+RGBA(:,1)-offsetvals(1) 2*offsetvals(2)*(handles.P(assocmask,1)/max(handles.P(assocmask,1)))+RGBA(:,2)-2*offsetvals(2) offsetvals(3)*(handles.P(assocmask,2)/max(handles.P(assocmask,2)))+RGBA(:,3)-offsetvals(3) ];

RGBA(:,1:3)=[2*offsetvals(1)*(handles.P(assocmask,2)/max(handles.P(assocmask,2)))+RGBA(:,1)-offsetvals(1)*2 3*offsetvals(2)*(handles.P(assocmask,1)/max(handles.P(assocmask,1)))+RGBA(:,2)-3*offsetvals(2) 3*offsetvals(3)*(handles.P(assocmask,2)/max(handles.P(assocmask,2)))+RGBA(:,3)-3*offsetvals(3) ];


% RGBA(:,4)=handles.origopac(i).*(1+handles.P(assocmask,2));
% %max(RGBA(:,4))
% RGBA(:,4)=(RGBA(:,4)/max(RGBA(:,4)))*handles.origopac(i);
% %max(RGBA(:,4))
 RGBA(:,1:3)=RGBA(:,1:3).*255;
% RGBA(RGBA>255)=255;


opval=handles.origopac(i).*(1+handles.P(assocmask,2));
if handles.selectallPressed && i==handles.totClusters-1
    opval=(opval/max(opval))*5*(handles.origopac(i)/max((handles.origopac)))+10+handles.alphaoffset(i);
else
    opval=(opval/max(opval))*120*(handles.origopac(i)/max((handles.origopac)))+80+handles.alphaoffset(i);
end
opval(opval>255)=255;
opval(opval<1)=1;
RGBA(:,4)=opval;




writeMat=horzcat(Ix,Iy,Iz,RGBA);


fid = fopen('../MyVTKData/select.bin', 'w','l');
fwrite(fid,1,'uint','l'); %select or cluster?
fwrite(fid, writeMat', 'float','l');
fclose(fid);
dummyval=1;
save('../MyVTKData/check.txt','dummyval','-ASCII');
while ~exist('../MyVTKData/doneren.txt','file')
    %pause(0.5);
end
delete('../MyVTKData/doneren.txt');
if exist('../MyVTKData/check.txt','file')
    delete('../MyVTKData/check.txt');
end

%delete('../MyVTKData/select.bin');

function hideorshowCluster(handles,i)


siz=size(handles.V);
assocmask=ismember(handles.assoc,find(handles.alreadySelected==i)); %assoc contains the association of each voxel to a node (1..642), find(ss) returns the index of who were selected, so this line returns logical 1 for voxels selected


selectInd=find(~handles.J(:));
selectInd=selectInd(assocmask);

[Ix,Iy,Iz]=ind2sub(siz,selectInd);
if ~handles.isclusteron(i)
    alphaval=ones(size(Ix,1),1);
else
    alphaval=handles.origopac(i).*(1+handles.P(assocmask,2));
    if  handles.selectallPressed && i==handles.totClusters-1
        
        alphaval=(alphaval/max(alphaval))*5*(handles.origopac(i)/max((handles.origopac)))+10+handles.alphaoffset(i);
        
    else
         alphaval=(alphaval/max(alphaval))*120*(handles.origopac(i)/max((handles.origopac)))+80+handles.alphaoffset(i);
    end
    alphaval(alphaval>255)=255;
    alphaval(alphaval<1)=1;
end

writeMat=horzcat(Ix,Iy,Iz,alphaval);
fid = fopen('../MyVTKData/select.bin', 'w','l');
fwrite(fid,0,'uint','l'); %select or cluster?
fwrite(fid,writeMat', 'float','l');
fclose(fid);
dummyval=1;
save('../MyVTKData/check.txt','dummyval','-ASCII');

while ~exist('../MyVTKData/doneren.txt','file')
    %pause(0.5);
end
delete('../MyVTKData/doneren.txt');
if exist('../MyVTKData/check.txt','file')
    delete('../MyVTKData/check.txt');
end

%delete('../MyVTKData/select.bin');


% --- Executes on selection change in huetemplatetype.
function huetemplatetype_Callback(hObject, eventdata, handles)

% hObject    handle to huetemplatetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns huetemplatetype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from huetemplatetype
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val};
    case 'X-Type'
        handles.hrange=[0.0 0.76]; % 76% of the hue ring
    case 'V-Type'
        handles.hrange=[0.0 0.26]; % 26% of the hue ring
    case 'T-Type'
        handles.hrange=[0.0 0.5]; % 50% of the hue ring
    case 'L-Type'
        handles.hrange=[0.0 0.385]; % 38.5% of the hue ring
        
end
handles.huevals=linspace(handles.hrange(1),handles.hrange(2),handles.totClusters-1);
[~,handles.huevals]=deal(fix(handles.huevals), handles.huevals-fix(handles.huevals));
%set(handles.hoffsetSlider, 'Max', 1.0-handles.hrange(2));
%set(handles.hoffsetSlider, 'Value', 0.0);
guidata(hObject, handles);
if handles.totClusters>1
    writeHSVforClusters(handles);
end
% --- Executes during object creation, after setting all properties.
function huetemplatetype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to huetemplatetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.


% --- Executes on slider movement.



% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
handles.isclusteron(1)=get(handles.checkbox1,'Value');
guidata(hObject, handles);

hideorshowCluster(handles,1);

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
handles.isclusteron(2)=get(handles.checkbox2,'Value');
guidata(hObject, handles);

hideorshowCluster(handles,2);
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
handles.isclusteron(3)=get(handles.checkbox3,'Value');
guidata(hObject, handles);
hideorshowCluster(handles,3);
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
handles.isclusteron(4)=get(handles.checkbox4,'Value');
guidata(hObject, handles);

hideorshowCluster(handles,4);
% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.isclusteron(5)=get(handles.checkbox5,'Value');
guidata(hObject, handles);

hideorshowCluster(handles,5);
% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in cplus1.
function cplus1_Callback(hObject, eventdata, handles)
% hObject    handle to cplus1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.alphaoffset(1)=handles.alphaoffset(1)+10;
guidata(hObject, handles);
hideorshowCluster(handles,1);
% --- Executes on button press in cminus1.
function cminus1_Callback(hObject, eventdata, handles)
% hObject    handle to cminus1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.alphaoffset(1)=handles.alphaoffset(1)-10;
guidata(hObject, handles);
hideorshowCluster(handles,1);

% --- Executes on button press in cplus2.
function cplus2_Callback(hObject, eventdata, handles)
% hObject    handle to cplus2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.alphaoffset(2)=handles.alphaoffset(2)+10;
guidata(hObject, handles);
hideorshowCluster(handles,2);
% --- Executes on button press in cminus2.
function cminus2_Callback(hObject, eventdata, handles)
% hObject    handle to cminus2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.alphaoffset(2)=handles.alphaoffset(2)-10;
guidata(hObject, handles);
hideorshowCluster(handles,2);

% --- Executes on button press in cplus3.
function cplus3_Callback(hObject, eventdata, handles)
% hObject    handle to cplus3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.alphaoffset(3)=handles.alphaoffset(3)+10;
guidata(hObject, handles);
hideorshowCluster(handles,3);


% --- Executes on button press in cminus3.
function cminus3_Callback(hObject, eventdata, handles)
% hObject    handle to cminus3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.alphaoffset(3)=handles.alphaoffset(3)-10;
guidata(hObject, handles);
hideorshowCluster(handles,3);

% --- Executes on button press in cplus4.
function cplus4_Callback(hObject, eventdata, handles)
% hObject    handle to cplus4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.alphaoffset(4)=handles.alphaoffset(4)+10;
guidata(hObject, handles);
hideorshowCluster(handles,4);

% --- Executes on button press in cminus4.
function cminus4_Callback(hObject, eventdata, handles)
% hObject    handle to cminus4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.alphaoffset(4)=handles.alphaoffset(4)-10;
guidata(hObject, handles);
hideorshowCluster(handles,4);

% --- Executes on button press in cplus5.
function cplus5_Callback(hObject, eventdata, handles)
% hObject    handle to cplus5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.alphaoffset(5)=handles.alphaoffset(5)+10;
guidata(hObject, handles);
hideorshowCluster(handles,5);

% --- Executes on button press in cminus5.
function cminus5_Callback(hObject, eventdata, handles)
% hObject    handle to cminus5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.alphaoffset(5)=handles.alphaoffset(5)-10;
guidata(hObject, handles);
hideorshowCluster(handles,5);


% --- Executes on button press in huerotate.
function huerotate_Callback(hObject, eventdata, handles)
% hObject    handle to huerotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.hrange(1)=handles.hrange(1)+0.05;
handles.hrange(2)=handles.hrange(2)+0.05;
handles.huevals=linspace(handles.hrange(1),handles.hrange(2),handles.totClusters-1);
[~,handles.huevals]=deal(fix(handles.huevals), handles.huevals-fix(handles.huevals));
guidata(hObject, handles);
writeHSVforClusters(handles);
