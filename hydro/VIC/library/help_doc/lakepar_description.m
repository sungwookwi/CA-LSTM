function varargout = lakepar_description(varargin)
% LAKEPAR_DESCRIPTION MATLAB code for lakepar_description.fig
%      LAKEPAR_DESCRIPTION, by itself, creates a new LAKEPAR_DESCRIPTION or raises the existing
%      singleton*.
%
%      H = LAKEPAR_DESCRIPTION returns the handle to a new LAKEPAR_DESCRIPTION or the handle to
%      the existing singleton*.
%
%      LAKEPAR_DESCRIPTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LAKEPAR_DESCRIPTION.M with the given input arguments.
%
%      LAKEPAR_DESCRIPTION('Property','Value',...) creates a new LAKEPAR_DESCRIPTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lakepar_description_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lakepar_description_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lakepar_description

% Last Modified by GUIDE v2.5 20-Oct-2015 10:39:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lakepar_description_OpeningFcn, ...
                   'gui_OutputFcn',  @lakepar_description_OutputFcn, ...
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


% --- Executes just before lakepar_description is made visible.
function lakepar_description_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lakepar_description (see VARARGIN)

% Choose default command line output for lakepar_description
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lakepar_description wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = lakepar_description_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
