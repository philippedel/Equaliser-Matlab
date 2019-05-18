function varargout = coursework_mat(varargin)
% COURSEWORK_MAT MATLAB code for coursework_mat.fig
%      COURSEWORK_MAT, by itself, creates a new COURSEWORK_MAT or raises the existing
%      singleton*.
%
%      H = COURSEWORK_MAT returns the handle to a new COURSEWORK_MAT or the handle to
%      the existing singleton*.
%
%      COURSEWORK_MAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COURSEWORK_MAT.M with the given input arguments.
%
%      COURSEWORK_MAT('Property','Value',...) creates a new COURSEWORK_MAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before coursework_mat_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to coursework_mat_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help coursework_mat

% Last Modified by GUIDE v2.5 12-Apr-2019 02:00:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @coursework_mat_OpeningFcn, ...
                   'gui_OutputFcn',  @coursework_mat_OutputFcn, ...
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


% --- Executes just before coursework_mat is made visible.
function coursework_mat_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to coursework_mat (see VARARGIN)

% Choose default command line output for coursework_mat
%handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
vol =  0.3;
set(handles.volume_slider,'value',vol);
handles.output = hObject;

guidata(hObject, handles);

% UIWAIT makes coursework_mat wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = coursework_mat_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in browse_button.
function browse_button_Callback(hObject, eventdata, handles)
% hObject    handle to browse_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname] = uigetfile({'*.mp3;*.wav;*.mp4'},'File Selector');
handles.fullpathname = strcat(pathname, filename);

set(handles.text_browse, 'String',handles.fullpathname) %showing fullpathname
guidata(hObject,handles)

% --- Executes on button press in play_button.
function play_button_Callback(hObject, eventdata, handles)
% hObject    handle to play_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global player;
play_equalizer(hObject, handles); 
play(player);
guidata(hObject,handles)

function play_equalizer(hObject, handles)
global player;
[handles.y,handles.Fs] = audioread(handles.fullpathname)
handles.Volume=get(handles.volume_slider,'value');

%NewStart = get(player, ' CurrentSample');
%handles.y = handles.y(NewStart :end,:);

handles.g1=get(handles.slider2,'value');
handles.g2=get(handles.slider4,'value');
handles.g3=get(handles.slider5,'value');
handles.g4=get(handles.slider6,'value');
handles.g5=get(handles.slider7,'value');
handles.g6=get(handles.slider8,'value');
handles.g7=get(handles.slider9,'value');
handles.g8=get(handles.slider10,'value');
handles.g9=get(handles.slider11,'value');
handles.g10=get(handles.slider12,'value');

%precision of 1 point
txt1 = sprintf('%.1f', handles.g1);
txt2 = sprintf('%.1f', handles.g2);
txt3 = sprintf('%.1f', handles.g3);
txt4 = sprintf('%.1f', handles.g4);
txt5 = sprintf('%.1f', handles.g5);
txt6 = sprintf('%.1f', handles.g6);
txt7 = sprintf('%.1f', handles.g7);
txt8 = sprintf('%.1f', handles.g8);
txt9 = sprintf('%.1f', handles.g9);
txt10 = sprintf('%.1f', handles.g10);

%setting the current value selected
set(handles.text2, 'String',txt1);
set(handles.text4, 'String',txt2);
set(handles.text5, 'String',txt3);
set(handles.text66, 'String',txt4);
set(handles.text7, 'String',txt5);
set(handles.text88, 'String',txt6);
set(handles.text9, 'String',txt7);
set(handles.text100, 'String',txt8);
set(handles.text111, 'String',txt9);
set(handles.text12, 'String',txt10);

handles.pop_menu = get(handles.popupmenu1,'Value');     %current filter value
handles.reverb_pop_menu = get(handles.reverb_popupmenu,'Value');   %current reverb value

handles.y_filter = handles.y;
switch handles.reverb_pop_menu
    case 1
    case 2
        %%%%%%%%
        %Impulse Response Reverberated Signal
        %%%%%%%%
        filename='impulse_demo.wav';
        [imp,Fsimp] = audioread(filename);
        handles.y_filter = fconv(handles.y_filter,imp);
        
        allItems = get(handles.popupmenu1,'string')     %getting the value of the filter
        selectedIndex = get(handles.popupmenu1,'Value')
        selectedItem = allItems{selectedIndex};
        
        title(['Impulse Response Reverberated Signal and ' selectedItem]);
    case 3
        %%%%%%%
        %moore reverb
        %%%%%%%
        % Call moorer reverb
        %set delay of each comb filter
        %set delay of each allpass filter in number of samples
        %Compute a random set of milliseconds and use sample rate
        rand('state',sum(100*clock))
        cd = floor(0.05*rand([1,6])*handles.Fs);

        % set gains of 6 comb pass filters
        g01 = 0.5*ones(1,6);
        %set feedback of each comb filter
        g02 = 0.5*ones(1,6);
        % set input cg and cg1 for moorer function see help moorer
        cg = g02./(1-g01);
        cg1 = g01;


        %set gain of allpass filter
        ag = 0.7;
        %set delay of allpass filter
        ad = 0.08*handles.Fs;
        %set direct signal gain
        k = 0.5;


        [handles.y_filter b a] = moorer(handles.y_filter,cg,cg1,cd,ag,ad,k);
        
        allItems = get(handles.popupmenu1,'string')     %getting the value of the filter
        selectedIndex = get(handles.popupmenu1,'Value')
        selectedItem = allItems{selectedIndex};
        
        title(['Moore Reverb and ' selectedItem]);
    case 4
        %%%%%%%
        %Schroeder reverb
        %%%%%%%

        % Call Schroeder1 reverb
        %set the number of allpass filters
        n = 6;
        %set the gain of the allpass filters
        g = 0.9;
        %set delay of each allpass filter in number of samples
        %Compute a random set of milliseconds and use sample rate
        rand('state',sum(100*clock))
        d = floor(0.05*rand([1,n])*handles.Fs);
        %set gain of direct signal
        k= 0.2;

        [handles.y_filter b a] = schroeder1(handles.y_filter,n,g,d,k);
        
        allItems = get(handles.popupmenu1,'string')     %getting the value of the filter
        selectedIndex = get(handles.popupmenu1,'Value')
        selectedItem = allItems{selectedIndex};
   
        title(['Schroeder Reverb and ' selectedItem]);
end

%check for wah wah

if (get(handles.radiobutton3,'Value') ~= 0)
    % damping factor
    % lower the damping factor the smaller the pass band
    damp = 0.05;
    
    % min and max centre cutoff frequency of variable bandpass filter
    minf=500;
    maxf=3000;
    
    % wah frequency, how many Hz per second are cycled through
    Fw = 2000;
   
    % change in centre frequency per sample (Hz)
    delta = Fw/handles.Fs;
    
    % create triangle wave of centre frequency values
    Fc=minf:delta:maxf;
    
    while(length(Fc) < length(handles.y_filter) )
        Fc= [ Fc (maxf:-delta:minf) ];
        Fc= [ Fc (minf:delta:maxf) ];
    end
    
    % trim tri wave to size of input
    Fc = Fc(1:length(handles.y_filter));
    
    % difference equation coefficients
    % must be recalculated each time Fc changes
    F1 = 2*sin((pi*Fc(1))/handles.Fs);
    
    % this dictates size of the pass bands
    Q1 = 2*damp;
    yh=zeros(size(handles.y_filter)); % create emptly out vectors
    yb=zeros(size(handles.y_filter));
    yl=zeros(size(handles.y_filter));
    
    % first sample, to avoid referencing of negative signals
    yh(1) = handles.y_filter(1);
    yb(1) = F1*yh(1);
    yl(1) = F1*yb(1);
    
    % apply difference equation to the sample
    for n=2:length(handles.y_filter),
        yh(n) = handles.y_filter(n) - yl(n-1) - Q1*yb(n-1);
        yb(n) = F1*yh(n) + yb(n-1);
        yl(n) = F1*yb(n) + yl(n-1);
        F1 = 2*sin((pi*Fc(n))/handles.Fs);
    end
    

    maxyb = max(abs(yb));
    handles.y_filter = yb/maxyb;
end

switch handles.pop_menu
    case 1
        %%%%%%
        %Base Shelving
        %%%%%%
        set(handles.text10, 'String','Hz');
        
        set(handles.text6, 'String','32');
        set(handles.text8, 'String','64');
        set(handles.text11, 'String','128');
        set(handles.text13, 'String','256');
        set(handles.text15, 'String','512');
        set(handles.text17, 'String','1200');
        set(handles.text19, 'String','2500');
        set(handles.text21, 'String','5500');
        set(handles.text23, 'String','11500');
        set(handles.text25, 'String','15000');
        G = 4;
        Q = 3;
        type = 'Base_Shelf';
        % %bandpass1
        fcb1 = 32;
        [b a] = shelving(handles.g1, fcb1, handles.Fs, Q, type);
        y1 = filter(b,a, handles.y_filter);
        % %bandpass2
        fcb2 = 64;
        [b a] = shelving(handles.g2, fcb2, handles.Fs, Q, type);
        y2 = filter(b,a, y1);
        % %bandpass3
        fcb3 = 128;
        [b a] = shelving(handles.g3, fcb3, handles.Fs, Q, type);
        y3 = filter(b,a, y2);
        % %bandpass4
        fcb4 = 256;
        [b a] = shelving(handles.g4, fcb4, handles.Fs, Q, type);
        y4 = filter(b,a, y3);
        % %bandpass5
        fcb5 = 512;
        [b a] = shelving(handles.g5, fcb5, handles.Fs, Q, type);
        y5 = filter(b,a, y4);
        % %bandpass6
        fcb6 = 1200;
        [b a] = shelving(handles.g6, fcb6, handles.Fs, Q, type);
        y6 = filter(b,a, y5);
        % %bandpass7
        fcb7 = 2500;
        [b a] = shelving(handles.g7, fcb7, handles.Fs, Q, type);
        y7 = filter(b,a, y6);
        % %bandpass8
        fcb8 = 5500;
        [b a] = shelving(handles.g8, fcb8, handles.Fs, Q, type);
        y8 = filter(b,a, y7);
        % %bandpass9
        fcb9 = 11500;
        [b a] = shelving(handles.g9, fcb9, handles.Fs, Q, type);
        y9 = filter(b,a, y7);
        % %bandpass10
        fcb10 = 15000;
        [b a] = shelving(handles.g10, fcb10, handles.Fs, Q, type);
        y10 = filter(b,a, y9);
        handles.y_filter = y10;
        
        cla
        hold on
        %plot(handles.y_filter,'b');
        plot(handles.y_filter,'b');
        plot(handles.y,'r');
        title('Bass Shelf Filter Equalised Signal');
        legend('BLUE is Filter','RED is Original');
    case 2
        %%%%%%
        %Treble Shelving
        %%%%%%
        set(handles.text10, 'String','Hz');
        
        set(handles.text6, 'String','32');
        set(handles.text8, 'String','64');
        set(handles.text11, 'String','128');
        set(handles.text13, 'String','256');
        set(handles.text15, 'String','512');
        set(handles.text17, 'String','1200');
        set(handles.text19, 'String','2500');
        set(handles.text21, 'String','5500');
        set(handles.text23, 'String','11500');
        set(handles.text25, 'String','15000');
        G = 4;
        Q = 3;
        type = 'Treble_Shelf';
        
        % %bandpass1
        fcb1 = 32;
        [b a] = shelving(handles.g1, fcb1, handles.Fs, Q, type);
        y1 = filter(b,a, handles.y_filter);
        % %bandpass2
        fcb2 = 64;
        [b a] = shelving(handles.g2, fcb2, handles.Fs, Q, type);
        y2 = filter(b,a, y1);
        % %bandpass3
        fcb3 = 128;
        [b a] = shelving(handles.g3, fcb3, handles.Fs, Q, type);
        y3 = filter(b,a, y2);
        % %bandpass4
        fcb4 = 256;
        [b a] = shelving(handles.g4, fcb4, handles.Fs, Q, type);
        y4 = filter(b,a, y3);
        % %bandpass5
        fcb5 = 512;
        [b a] = shelving(handles.g5, fcb5, handles.Fs, Q, type);
        y5 = filter(b,a, y4);
        % %bandpass6
        fcb6 = 1200;
        [b a] = shelving(handles.g6, fcb6, handles.Fs, Q, type);
        y6 = filter(b,a, y5);
        % %bandpass7
        fcb7 = 2500;
        [b a] = shelving(handles.g7, fcb7, handles.Fs, Q, type);
        y7 = filter(b,a, y6);
        % %bandpass8
        fcb8 = 5500;
        [b a] = shelving(handles.g8, fcb8, handles.Fs, Q, type);
        y8 = filter(b,a, y7);
        % %bandpass9
        fcb9 = 11500;
        [b a] = shelving(handles.g9, fcb9, handles.Fs, Q, type);
        y9 = filter(b,a, y7);
        % %bandpass10
        fcb10 = 15000;
        [b a] = shelving(handles.g10, fcb10, handles.Fs, Q, type);
        y10 = filter(b,a, y9);
        handles.y_filter = y10;
        cla
        hold on
        plot(handles.y_filter,'b');
        plot(handles.y,'r');
        title('Treble Shelf Filter Equalised Signal');
        legend('BLUE is Filter','RED is Original');
    case 3
        set(handles.text10, 'String','');
        
        set(handles.text6, 'String','0-200 hz');
        set(handles.text8, 'String','201-400 hz');
        set(handles.text11, 'String','401-800 hz');
        set(handles.text13, 'String','801-1500 hz');
        set(handles.text15, 'String','1.5-3 Khz');
        set(handles.text17, 'String','3-5 Khz');
        set(handles.text19, 'String','5-7 Khz');
        set(handles.text21, 'String','7-10 Khz');
        set(handles.text23, 'String','10-15 Khz');
        set(handles.text25, 'String','>15 Khz');
        
        cut_off=200; %cut off low pass dalama Hz
        orde=16;
        a=fir1(orde,cut_off/(handles.Fs/2),'low');
        %lowpass, bandpass, or multiband FIR filter with linear phase. The filter type depends on the number of elements of Wn .
        y1=handles.g1/3*filter(a,1,handles.y_filter);

        % %bandpass1
        f1=201; 
        f2=400; 
        b1=fir1(orde,[f1/(handles.Fs/2) f2/(handles.Fs/2)],'bandpass');
        y2=handles.g2/3*filter(b1,1,handles.y);
        % 
        % %bandpass2
        f3=401;
        f4=800;
        b2=fir1(orde,[f3/(handles.Fs/2) f4/(handles.Fs/2)],'bandpass');
        y3=handles.g3/3*filter(b2,1,handles.y);
        % 
        % %bandpass3
         f4=801;
        f5=1500;
         b3=fir1(orde,[f4/(handles.Fs/2) f5/(handles.Fs/2)],'bandpass');
         y4=handles.g4/3*filter(b3,1,handles.y);
        % 
        % %bandpass4
         f5=1501;
        f6=3000;
         b4=fir1(orde,[f5/(handles.Fs/2) f6/(handles.Fs/2)],'bandpass');
         y5=handles.g5/3*filter(b4,1,handles.y);
        % 
        % %bandpass5
          f7=3001;
        f8=5000;
          b5=fir1(orde,[f7/(handles.Fs/2) f8/(handles.Fs/2)],'bandpass');
          y6=handles.g6/3*filter(b5,1,handles.y);
        % 
        % %bandpass6
          f9=5001;
        f10=7000;
          b6=fir1(orde,[f9/(handles.Fs/2) f10/(handles.Fs/2)],'bandpass');
          y7=handles.g7/3*filter(b6,1,handles.y);
        % 
        % %bandpass7
          f11=7001;
        f12=10000;
          b7=fir1(orde,[f11/(handles.Fs/2) f12/(handles.Fs/2)],'bandpass');
          y8=handles.g8/3*filter(b7,1,handles.y);
        % 
         % %bandpass8
          f13=10001;
        f14=15000;
          b8=fir1(orde,[f13/(handles.Fs/2) f14/(handles.Fs/2)],'bandpass');
          y9=handles.g9/3*filter(b8,1,handles.y);
        % 
         %highpass
        cut_off2=15000;
        c=fir1(orde,cut_off2/(handles.Fs/2),'high');
        y10=handles.g10/3*filter(c,1,handles.y);
        
        handles.y_filter = y1+y2+y3+y4+y5+y6+y7+y8+y9+y10;
        cla
        hold on
        plot(handles.y_filter,'b')
        plot(handles.y,'r');
        title('FIR Filter Equalised Signal');
        legend('BLUE is Filter','RED is Original');
        
        otherwise
end

%hold on
%ten = 0:10;
%ten_slider = [handles.g1 handles.g2 handles.g3 handles.g4 handles.g6 handles.g7 handles.g8 handles.g9 handles.g10];

%plot(ten,ten_slider,'g');

%handles.y_filter = y1+y2+y3+y4+y5+y6+y7+y8+y9+y10;
player = audioplayer(handles.Volume/10*handles.y_filter, handles.Fs);
guidata(hObject,handles)

function [b, a]  = shelving(G, fc, fs, Q, type)
%
% Derive coefficients for a shelving filter with a given amplitude and
% cutoff frequency.  All coefficients are calculated as described in 
% Zolzer's DAFX book (p. 50 -55).  
%
% Usage:     [B,A] = shelving(G, Fc, Fs, Q, type);
%
%            G is the logrithmic gain (in dB)
%            FC is the center frequency
%            Fs is the sampling rate
%            Q adjusts the slope be replacing the sqrt(2) term
%            type is a character string defining filter type
%                 Choices are: 'Base_Shelf' or 'Treble_Shelf'
%
% Author:    Jeff Tackett 08/22/05
%

%Error Check
if((strcmp(type,'Base_Shelf') ~= 1) && (strcmp(type,'Treble_Shelf') ~= 1))
    error(['Unsupported Filter Type: ' type]);
end

K = tan((pi * fc)/fs);
V0 = 10^(G/20);
root2 = 1/Q; %sqrt(2)

%Invert gain if a cut
if(V0 < 1)
    V0 = 1/V0;
end

%%%%%%%%%%%%%%%%%%%%
%    BASE BOOST
%%%%%%%%%%%%%%%%%%%%
if(( G > 0 ) & (strcmp(type,'Base_Shelf')))
   
    b0 = (1 + sqrt(V0)*root2*K + V0*K^2) / (1 + root2*K + K^2);
    b1 =             (2 * (V0*K^2 - 1) ) / (1 + root2*K + K^2);
    b2 = (1 - sqrt(V0)*root2*K + V0*K^2) / (1 + root2*K + K^2);
    a1 =                (2 * (K^2 - 1) ) / (1 + root2*K + K^2);
    a2 =             (1 - root2*K + K^2) / (1 + root2*K + K^2);

%%%%%%%%%%%%%%%%%%%%
%    BASE CUT
%%%%%%%%%%%%%%%%%%%%
elseif (( G < 0 ) & (strcmp(type,'Base_Shelf')))
    
    b0 =             (1 + root2*K + K^2) / (1 + root2*sqrt(V0)*K + V0*K^2);
    b1 =                (2 * (K^2 - 1) ) / (1 + root2*sqrt(V0)*K + V0*K^2);
    b2 =             (1 - root2*K + K^2) / (1 + root2*sqrt(V0)*K + V0*K^2);
    a1 =             (2 * (V0*K^2 - 1) ) / (1 + root2*sqrt(V0)*K + V0*K^2);
    a2 = (1 - root2*sqrt(V0)*K + V0*K^2) / (1 + root2*sqrt(V0)*K + V0*K^2);

%%%%%%%%%%%%%%%%%%%%
%   TREBLE BOOST
%%%%%%%%%%%%%%%%%%%%
elseif (( G > 0 ) & (strcmp(type,'Treble_Shelf')))

    b0 = (V0 + root2*sqrt(V0)*K + K^2) / (1 + root2*K + K^2);
    b1 =             (2 * (K^2 - V0) ) / (1 + root2*K + K^2);
    b2 = (V0 - root2*sqrt(V0)*K + K^2) / (1 + root2*K + K^2);
    a1 =              (2 * (K^2 - 1) ) / (1 + root2*K + K^2);
    a2 =           (1 - root2*K + K^2) / (1 + root2*K + K^2);

%%%%%%%%%%%%%%%%%%%%
%   TREBLE CUT
%%%%%%%%%%%%%%%%%%%%

elseif (( G < 0 ) & (strcmp(type,'Treble_Shelf')))

    b0 =               (1 + root2*K + K^2) / (V0 + root2*sqrt(V0)*K + K^2);
    b1 =                  (2 * (K^2 - 1) ) / (V0 + root2*sqrt(V0)*K + K^2);
    b2 =               (1 - root2*K + K^2) / (V0 + root2*sqrt(V0)*K + K^2);
    a1 =             (2 * ((K^2)/V0 - 1) ) / (1 + root2/sqrt(V0)*K + (K^2)/V0);
    a2 = (1 - root2/sqrt(V0)*K + (K^2)/V0) / (1 + root2/sqrt(V0)*K + (K^2)/V0);

%%%%%%%%%%%%%%%%%%%%
%   All-Pass
%%%%%%%%%%%%%%%%%%%%
else
    b0 = V0;
    b1 = 0;
    b2 = 0;
    a1 = 0;
    a2 = 0;
end

%return values
a = [  1, a1, a2];
b = [ b0, b1, b2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%until here in file
function [y,b,a]=schroeder1(x,n,g,d,k)
%This is a reverberator based on Schroeder's design which consists of n all pass filters in series.
%
%The structure is:  [y,b,a] = schroeder1(x,n,g,d,k) 
%
%where x = the input signal
%      n = the number of allpass filters	
%      g = the gain of the allpass filters (this should be less than 1 for stability)
%      d = a vector which contains the delay length of each allpass filter
%      k = the gain factor of the direct signal
%      y = the output signal
%      b = the numerator coefficients of the transfer function
%      a = the denominator coefficients of the transfer function
%
% note: Make sure that d is the same length as n.
%
%
% Gautham J. Mysore - gauthamjm@yahoo.com
%



% send the input signal through the first allpass filter
[y,b,a] = allpass(x,g,d(1));

% send the output of each allpass filter to the input of the next allpass filter
for i = 2:n,
   [y,b1,a1] = allpass(y,g,d(i));
   [b,a] = seriescoefficients(b1,a1,b,a);
end   

% add the scaled direct signal
y = y + k*x;

% normalize the output signal
y = y/max(y);

function [y]=fconv(x, h)
%FCONV Fast Convolution
%   [y] = FCONV(x, h) convolves x and h, and normalizes the output  
%         to +-1.
%
%      x = input vector
%      h = input vector
% 
%      See also CONV
%
%   NOTES:
%
%   1) I have a short article explaining what a convolution is.  It
%      is available at http://stevem.us/fconv.html.
%
%
%Version 1.0
%Coded by: Stephen G. McGovern, 2003-2004.

Ly=length(x)+length(h)-1;  % 
Ly2=pow2(nextpow2(Ly));    % Find smallest power of 2 that is > Ly
X=fft(x, Ly2);		   % Fast Fourier transform
H=fft(h, Ly2);	           % Fast Fourier transform
Y=X.*H;        	           % 
y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly);               % Take just the first N elements
y=y/max(abs(y));           % Normalize the output

function [y,b,a]=moorer(x,cg,cg1,cd,ag,ad,k)
%This is a reverberator based on Moorer's design which consists of 6 parallel feedback comb filters 
%(each with a low pass filter in the feedback loop) in series with an all pass filter.
%
%The structure is:  [y,b,a] = moorer(x,cg,cg1,cd,ag,ad,k)
%
%where x = the input signal
%      cg = a vector of length 6 which contains g2/(1-g1) (this should be less than 1 for stability),
%           where g2 is the feedback gain of each of the comb filters and g1 is from the following parameter 
%      cg1 = a vector of length 6 which contains the gain of the low pass filters in the feedback loop of
%            each of the comb filters (should be less than 1 for stability)
%      cd = a vector of length 6 which contains the delay of each of the comb filters 
%      ag = the gain of the allpass filter (this should be less than 1 for stability)
%      ad = the delay of the allpass filter 
%      k = the gain factor of the direct signal
%      y = the output signal
%      b = the numerator coefficients of the transfer function
%      a = the denominator coefficients of the transfer function
%
%
% Gautham J. Mysore - gauthamjm@yahoo.com
%


% send the input to each of the 6 comb filters separately
[outcomb1,b1,a1] = lpcomb(x,cg(1),cg1(1),cd(1));
[outcomb2,b2,a2] = lpcomb(x,cg(2),cg1(2),cd(2));
[outcomb3,b3,a3] = lpcomb(x,cg(3),cg1(3),cd(3));
[outcomb4,b4,a4] = lpcomb(x,cg(4),cg1(4),cd(4));
[outcomb5,b5,a5] = lpcomb(x,cg(5),cg1(5),cd(5));
[outcomb6,b6,a6] = lpcomb(x,cg(6),cg1(6),cd(6));

% sum the ouptut of the 6 comb filters
apinput = outcomb1 + outcomb2 + outcomb3 + outcomb4 + outcomb5 + outcomb6; 

%find the combined filter coefficients of the the comb filters
[b,a]=parallelcoefficients(b1,a1,b2,a2);
[b,a]=parallelcoefficients(b,a,b3,a3);
[b,a]=parallelcoefficients(b,a,b4,a4);
[b,a]=parallelcoefficients(b,a,b5,a5);
[b,a]=parallelcoefficients(b,a,b6,a6);

% send the output of the comb filters to the allpass filter
[y,b7,a7] = allpass(apinput,ag,ad);

%find the combined filter coefficients of the the comb filters in series with the allpass filters
[b,a]=seriescoefficients(b,a,b7,a7);

% add the scaled direct signal
y = y + k*x;

% normalize the output signal
y = y/max(y);

% --- Executes on button press in pause_button.
function pause_button_Callback(hObject, eventdata, handles)
% hObject    handle to pause_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global player;
pause(player);
guidata(hObject,handles)

function [y,b,a]=lpcomb(x,g,g1,d)

%This is a feedback comb filter with a low pass filter in the feedback.
%
%The structure is:  [y,b,a] = lpcomb(x,g,g1,d)
%
%where x = the input signal
%      g = g2/(1-g1), where g2 is the feedback gain of the comb filter (this should be less than 1 for stability)
%      g1 = the feedback gain of the low pass filter (this should be less than 1 for stability)
%      d = the delay length
%      y = the output signal
%      b = the numerator coefficients of the transfer function
%      a = the denominator coefficients of the transfer function
%
%
% Gautham J. Mysore - gauthamjm@yahoo.com
%


%If g is more than 1, set it to 0.7 .
if g>=1
   g=0.7;
end   

%If the low pass feedback gain is more than 1, set it to 0.7 .
if g>=1
   g=0.7;
end   

%Set the b and a coefficients of the transfer function depending on g, g1 and d.
b=[zeros(1,d) 1 -g1];
a=[1 -g1 zeros(1,d-2) -g*(1-g1)];

%filter the input signal 
y=filter(b,a,x);

function [b,a]=parallelcoefficients(b1,a1,b2,a2)

%This function gives the filter coefficients of the parallel connection of two filters.
%
%The structure is:  [b,a] = parallelcoefficients(b1,a1,b2,a2)
%
%where b1 = the numerator coefficients of the 1st transfer function
%      a1 = the denominator coefficients of the 1st transfer function
%      b2 = the numerator coefficients of the 2nd transfer function
%      a2 = the denominator coefficients of the 2nd transfer function
%      b = the numerator coefficients of the composite transfer function
%      a = the denominator coefficients of the composite transfer function
%
%
% Gautham J. Mysore - gauthamjm@yahoo.com
%

b = conv(b1,a2) + conv(b2,a1);

a = conv(a1,a2);

function [y,b,a]=allpass(x,g,d)

%This is an allpass filter function.
%
%The structure is:  [y,b,a] = allpass(x,g,d)
%
%where x = the input signal
%      g = the feedforward gain (the feedback gain is the negative of this) (this should be less than 1 for stability)
%      d = the delay length
%      y = the output signal
%      b = the numerator coefficients of the transfer function
%      a = the denominator coefficients of the transfer function
%
%
% Gautham J. Mysore - gauthamjm@yahoo.com
%


%If the feedback gain is more than 1, set it to 0.7 .
if g>=1
   g=0.7;
end   

%Set the b and a coefficients of the transfer function depending on g and d.
b=[g zeros(1,d-1) 1];
a=[1 zeros(1,d-1) g];

%filter the input signal 
y=filter(b,a,x);

function [b,a]=seriescoefficients(b1,a1,b2,a2)

%This function gives the filter coefficients of the series connection of two filters.
%
%The structure is:  [b,a] = seriescoefficients(b1,a1,b2,a2)
%
%where b1 = the numerator coefficients of the 1st transfer function
%      a1 = the denominator coefficients of the 1st transfer function
%      b2 = the numerator coefficients of the 2nd transfer function
%      a2 = the denominator coefficients of the 2nd transfer function
%      b = the numerator coefficients of the composite transfer function
%      a = the denominator coefficients of the composite transfer function
%
%
% Gautham J. Mysore - gauthamjm@yahoo.com
%

b = conv(b1,b2);

a = conv(a1,a2);





% --- Executes on button press in resume_button.
function resume_button_Callback(hObject, eventdata, handles)
% hObject    handle to resume_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global player;
resume(player);
guidata(hObject,handles)


% --- Executes on button press in stop_button.
function stop_button_Callback(hObject, eventdata, handles)
% hObject    handle to stop_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global player;
play_equalizer(hObject, handles); 
stop(player);
guidata(hObject,handles)


% --- Executes on slider movement.
function volume_slider_Callback(hObject, eventdata, handles)
% hObject    handle to volume_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global player;
play_equalizer(hObject, handles); 
play(player);



% --- Executes during object creation, after setting all properties.
function volume_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to volume_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider8_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider9_Callback(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global player;
play_equalizer(hObject, handles); 
play(player);

% --- Executes during object creation, after setting all properties.
function slider9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider10_Callback(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes during object creation, after setting all properties.
function slider10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider11_Callback(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes during object creation, after setting all properties.
function slider11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider12_Callback(hObject, eventdata, handles)
% hObject    handle to slider12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global player;
play_equalizer(hObject, handles); 
play(player);

% --- Executes during object creation, after setting all properties.
function slider12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
handles = guidata(hObject);  % Update!
handles.pop_menu = get(handles.popupmenu1,'Value');
guidata(hObject, handles);  % Update!
global player;
play_equalizer(hObject, handles); 
play(player);





% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Rock_button.
function Rock_button_Callback(hObject, eventdata, handles)
% hObject    handle to Rock_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g1 = -7.5;
g2 = -8.6;
g3 = -3.6;
g4 = -2.7;
g5 =  2.1;
g6 = 5;
g7 = 5.5;
g8= 4.8;
g9 =4.8;
g10 = 5.1;
set(handles.slider2,'value',g1);
set(handles.slider4,'value',g2);
set(handles.slider5,'value',g3);
set(handles.slider6,'value',g4);
set(handles.slider7,'value',g5);
set(handles.slider8,'value',g6);
set(handles.slider9,'value',g7);
set(handles.slider10,'value',g8);
set(handles.slider11,'value',g9);
set(handles.slider12,'value',g10);

set(handles.text2, 'String',g1);
set(handles.text4, 'String',g2);
set(handles.text5, 'String',g3);
set(handles.text66, 'String',g4);
set(handles.text7, 'String',g5);
set(handles.text88, 'String',g6);
set(handles.text9, 'String',g7);
set(handles.text100, 'String',g8);
set(handles.text111, 'String',g9);
set(handles.text12, 'String',g10);

global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes on button press in Classical_button.
function Classical_button_Callback(hObject, eventdata, handles)
% hObject    handle to Classical_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g1 = 0;
g2 = 0;
g3 = 0;
g4 = 0;
g5 =  2;
g6 = 1;
g7 = -0.3;
g8= -5.7;
g9 = -6;
g10 = -8.1;

set(handles.slider2,'value',g1);
set(handles.slider4,'value',g2);
set(handles.slider5,'value',g3);
set(handles.slider6,'value',g4);
set(handles.slider7,'value',g5);
set(handles.slider8,'value',g6);
set(handles.slider9,'value',g7);
set(handles.slider10,'value',g8);
set(handles.slider11,'value',g9);
set(handles.slider12,'value',g10);

set(handles.text2, 'String',g1);
set(handles.text4, 'String',g2);
set(handles.text5, 'String',g3);
set(handles.text66, 'String',g4);
set(handles.text7, 'String',g5);
set(handles.text88, 'String',g6);
set(handles.text9, 'String',g7);
set(handles.text100, 'String',g8);
set(handles.text111, 'String',g9);
set(handles.text12, 'String',g10);

global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes on button press in Pop_button.
function Pop_button_Callback(hObject, eventdata, handles)
% hObject    handle to Pop_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g1 = -1.5;
g2 = 3.9;
g3 = 5.4;
g4 = 4.5;
g5 =  0.9;
g6 = -1.5;
g7 = -1.8;
g8= -2.1;
g9 = -2.1;
g10 = -0.3;

set(handles.slider2,'value',g1);
set(handles.slider4,'value',g2);
set(handles.slider5,'value',g3);
set(handles.slider6,'value',g4);
set(handles.slider7,'value',g5);
set(handles.slider8,'value',g6);
set(handles.slider9,'value',g7);
set(handles.slider10,'value',g8);
set(handles.slider11,'value',g9);
set(handles.slider12,'value',g10);

set(handles.text2, 'String',g1);
set(handles.text4, 'String',g2);
set(handles.text5, 'String',g3);
set(handles.text66, 'String',g4);
set(handles.text7, 'String',g5);
set(handles.text88, 'String',g6);
set(handles.text9, 'String',g7);
set(handles.text100, 'String',g8);
set(handles.text111, 'String',g9);
set(handles.text12, 'String',g10);

global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes on button press in Techno_button.
function Techno_button_Callback(hObject, eventdata, handles)
% hObject    handle to Techno_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g1 = 2.8;
g2 = 2.2;
g3 = 0.5;
g4 = -2.4;
g5 =  -3.3;
g6 = -1.5;
g7 = 1.5;
g8= 4.1;
g9 = 4.7;
g10 = 4.4;

set(handles.slider2,'value',g1);
set(handles.slider4,'value',g2);
set(handles.slider5,'value',g3);
set(handles.slider6,'value',g4);
set(handles.slider7,'value',g5);
set(handles.slider8,'value',g6);
set(handles.slider9,'value',g7);
set(handles.slider10,'value',g8);
set(handles.slider11,'value',g9);
set(handles.slider12,'value',g10);

set(handles.text2, 'String',g1);
set(handles.text4, 'String',g2);
set(handles.text5, 'String',g3);
set(handles.text66, 'String',g4);
set(handles.text7, 'String',g5);
set(handles.text88, 'String',g6);
set(handles.text9, 'String',g7);
set(handles.text100, 'String',g8);
set(handles.text111, 'String',g9);
set(handles.text12, 'String',g10);

global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes on button press in Blues_button.
function Blues_button_Callback(hObject, eventdata, handles)
% hObject    handle to Blues_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g1 = 1;
g2 = 0.5;
g3 = -0.6;
g4 = -2;
g5 =  -1.4;
g6 = -1;
g7 = -1;
g8= 2;
g9 = 3;
g10 = 2;

set(handles.slider2,'value',g1);
set(handles.slider4,'value',g2);
set(handles.slider5,'value',g3);
set(handles.slider6,'value',g4);
set(handles.slider7,'value',g5);
set(handles.slider8,'value',g6);
set(handles.slider9,'value',g7);
set(handles.slider10,'value',g8);
set(handles.slider11,'value',g9);
set(handles.slider12,'value',g10);

set(handles.text2, 'String',g1);
set(handles.text4, 'String',g2);
set(handles.text5, 'String',g3);
set(handles.text66, 'String',g4);
set(handles.text7, 'String',g5);
set(handles.text88, 'String',g6);
set(handles.text9, 'String',g7);
set(handles.text100, 'String',g8);
set(handles.text111, 'String',g9);
set(handles.text12, 'String',g10);

global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes on button press in Reset_button.
function Reset_button_Callback(hObject, eventdata, handles)
% hObject    handle to Reset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g1 = 0;
g2 = 0;
g3 = 0;
g4 = 0;
g5 =  0;
g6 = 0;
g7 = 0;
g8= 0;
g9 = 0;
g10 = 0;

set(handles.slider2,'value',g1);
set(handles.slider4,'value',g2);
set(handles.slider5,'value',g3);
set(handles.slider6,'value',g4);
set(handles.slider7,'value',g5);
set(handles.slider8,'value',g6);
set(handles.slider9,'value',g7);
set(handles.slider10,'value',g8);
set(handles.slider11,'value',g9);
set(handles.slider12,'value',g10);

set(handles.text2, 'String',g1);
set(handles.text4, 'String',g2);
set(handles.text5, 'String',g3);
set(handles.text66, 'String',g4);
set(handles.text7, 'String',g5);
set(handles.text88, 'String',g6);
set(handles.text9, 'String',g7);
set(handles.text100, 'String',g8);
set(handles.text111, 'String',g9);
set(handles.text12, 'String',g10);

global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes on button press in Costum_button.
function Costum_button_Callback(hObject, eventdata, handles)
% hObject    handle to Costum_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%open file prevoisly created with save button
fileID = fopen('saved.txt','r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
fclose(fileID);

set(handles.slider2,'value',A(1));
set(handles.slider4,'value',A(2));
set(handles.slider5,'value',A(3));
set(handles.slider6,'value',A(4));
set(handles.slider7,'value',A(5));
set(handles.slider8,'value',A(6));
set(handles.slider9,'value',A(7));
set(handles.slider10,'value',A(8));
set(handles.slider11,'value',A(9));
set(handles.slider12,'value',A(10));

set(handles.text2, 'String',A(1));
set(handles.text4, 'String',A(2));
set(handles.text5, 'String',A(3));
set(handles.text66, 'String',A(4));
set(handles.text7, 'String',A(5));
set(handles.text88, 'String',A(6));
set(handles.text9, 'String',A(7));
set(handles.text100, 'String',A(8));
set(handles.text111, 'String',A(9));
set(handles.text12, 'String',A(10));



global player;
play_equalizer(hObject, handles); 
play(player);



% --- Executes on button press in Save_costum_button.
function Save_costum_button_Callback(hObject, eventdata, handles)
% hObject    handle to Save_costum_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.g11=get(handles.slider2,'value');
handles.g22=get(handles.slider4,'value');
handles.g33=get(handles.slider5,'value');
handles.g44=get(handles.slider6,'value');
handles.g55=get(handles.slider7,'value');
handles.g66=get(handles.slider8,'value');
handles.g77=get(handles.slider9,'value');
handles.g88=get(handles.slider10,'value');
handles.g99=get(handles.slider11,'value');
handles.g100=get(handles.slider12,'value');

%save data in txt file for later use, if program is closed
fileID = fopen('saved.txt','w'); %open file
fprintf(fileID,'%.1f %f\n',handles.g11);
fprintf(fileID,'%.1f %f\n',handles.g22);
fprintf(fileID,'%.1f %f\n',handles.g33);
fprintf(fileID,'%.1f %f\n',handles.g44);
fprintf(fileID,'%.1f %f\n',handles.g55);
fprintf(fileID,'%.1f %f\n',handles.g66);
fprintf(fileID,'%.1f %f\n',handles.g77);
fprintf(fileID,'%.1f %f\n',handles.g88);
fprintf(fileID,'%.1f %f\n',handles.g99);
fprintf(fileID,'%.1f %f\n',handles.g100);

fclose(fileID); %close file

Costum_button_Callback(hObject, eventdata, handles);

guidata(hObject, handles);


% --- Executes on button press in Save_fig_button.
function Save_fig_button_Callback(hObject, eventdata, handles)
% hObject    handle to Save_fig_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in spectro_button.
function spectro_button_Callback(hObject, eventdata, handles)
% hObject    handle to spectro_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mono = (handles.y_filter + handles.Fs)/2;
figure(1);
spectrogram(mono,1024,120,2048,'power','yaxis');
set(gca,'YScale','log');
colormap('winter');
xlabel('Time')
ylabel('Frequency (Log Scale)')

shg;


% --- Executes on button press in Record_button.
function Record_button_Callback(hObject, eventdata, handles)
% hObject    handle to Record_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
audiowrite('saved_song.wav', handles.y_filter, handles.Fs);


% --- Executes on selection change in reverb_popupmenu.
function reverb_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to reverb_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns reverb_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from reverb_popupmenu
global player;
play_equalizer(hObject, handles); 
play(player);


% --- Executes during object creation, after setting all properties.
function reverb_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reverb_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
global player;
play_equalizer(hObject, handles); 
play(player);