function varargout = pivot65(varargin)
% PIVOT65 Application M-file for pivot65.fig
% type pivot65 at the command line, making sure
% the active directory contains PIVOT's files

% Edit the above text to modify the response to help pivot65

% Last Modified by GUIDE v2.5 08-May-2006 21:47:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',          mfilename, ...
    'gui_Singleton',     gui_Singleton, ...
    'gui_OpeningFcn',    @pivot65_OpeningFcn, ...
    'gui_OutputFcn',     @pivot65_OutputFcn, ...
    'gui_LayoutFcn',     [], ...
    'gui_Callback',      []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    varargout{1:nargout} = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before pivot65 is made visible.
function pivot65_OpeningFcn(hObject, eventdata, handles, varargin, nargin)
global fileSpecified
global nargin
if nargin
    fileName = varargin{1};
    fileParams = readInputFile(fileName);
    fileSpecified = 1;
    fieldFiller;
end

global module
if fileSpecified
    set(handles,'Value','')
    ...
end

set(handles.edit9,'Visible','Off');
set(handles.pushbutton4,'Visible','Off');
set(handles.edit13,'Visible','Off');
set(handles.edit9,'String','');
set(handles.slider2,'Visible','Off');
set(handles.xSlider,'Visible','Off');
set(handles.xField,'Visible','Off');
set(handles.tSlider,'Visible','Off');
set(handles.tField,'Visible','Off');


% Choose default command line output for pivot65
handles.output = {};

% Update handles structure
guidata(hObject, handles);

% Selects axes and generates default plots for them
axes(handles.surfacePlots);
surf(zeros(100,100));
xlabel('position'),ylabel('time'),zlabel('temperature');
axes(handles.linePlots);
plot(zeros(1,10));
grid on;

module = '';
eqWindow = eqInfo65
movegui(eqWindow, 'northwest');

% --- Outputs from this function are returned to the command line.
function varargout = pivot65_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% --------------------------------------------------------------------
function varargout = plot_button_Callback(h, eventdata, handles, varargin)
syms x; syms t; syms y;

global module u x_coords t_coords y_coords xt tx audioON numModes ICF timeSpan pivotMovie L uu uu_Big xx yy tt;
global A_coefficients B_coefficients A_basis B_basis userIC_u userIC_dudt upperSelection;
global T X_U tt;
global xLeft xRight yLeft yRight;
clear global pivotMovie;

set(handles.pushbutton4,'Visible','Off');
evalBar = waitbar(1,'Evaluating.  Please wait.');
% Get user input from GUI
set(handles.upperPlot,'Value',2);
upperSelection = 'Surface Plot';

numModes = str2double(get(handles.numModes,'String'));
k= str2double(get(handles.f2_input,'String'));
L = str2double(get(handles.L_input,'String'));
H = str2double(get(handles.Height,'String'));
alpha = str2double(get(handles.alpha_edit,'String'));
beta = str2double(get(handles.beta_edit,'String'));
b = str2double(get(handles.b_input,'String'));

IC_u = sym(get(handles.IC,'String'));
if strcmp(char(IC_u),'USER-DEFINED')
    IC_u = userIC_u;
end
timeSpan = str2double(get(handles.tSpan,'String'));

if strcmp(get(handles.IC_derivative,'String'),'USER-DEFINED')
    IC_dudt = userIC_dudt;
else
    try
        IC_dudt = sym(get(handles.IC_derivative,'String'));
    catch
        IC_dudt = '0'
    end
end


if or(module == '1DHEAT',module == '1DChar')
    try
        source = sym(get(handles.sourceFxn,'String'));
    catch
        source = str2num(get(handles.sourceFxn,'String'));
    end
else
    source = str2num(get(handles.sourceFxn,'String'))
    if source ~= 0
        source = 0;
        warndlg('Source terms are only currently available for the 1-D Heat Equation. The source has been set to 0');
    end
end

if or(strcmp(get(handles.alpha_edit,'String'),'Edited'),source ~= 0)
    if strcmp(get(handles.alpha_edit,'String'),'Edited')
        try
            leftX = xLeft(3)/xLeft(1);
        end

        try
            rightX = xRight(3)/xRight(1);
        end
        PDEQ = 'NH';
    else
        PDEQ = 'NH';

        leftX = 0;
        rightX = 0;
    end
else
    PDEQ = 'HG';
end

if module == '1DHEAT'
    if PDEQ == 'HG'
        [u, A_coefficients, B_coefficients, A_basis, B_basis] = HeatEquation(numModes, k, L, alpha, beta, IC_u);
        x_coords = linspace(0,L); t_coords = linspace(0,timeSpan);
        [xt,tx] = meshgrid(x_coords, t_coords);

        if diff(u,x) ~= 0
            uu = subs(u,{x,t}, {xt,tx});
        else
            uu = double(u) * ones(size(xt));
        end

        uu_Big = max(max(abs(uu)));
        axes(handles.surfacePlots)
        surf(xt,tx,uu);
        xlabel('position'); ylabel('time'); zlabel('temperature');
        set(handles.upperPlot,'Value',2);
    else
        try
            [u, A_coefficients, B_coefficients, A_basis, B_basis, v, u_Equilibrium] = nonHomogHeat(numModes,k,L, alpha, beta, IC_u, source, leftX, rightX);
            x_coords = linspace(0,L); t_coords = linspace(0,timeSpan);
            [xt,tx] = meshgrid(x_coords, t_coords);

            if diff(u,x) ~= 0
                uu = subs(u,{x,t}, {xt,tx});
            else
                uu = double(u) * ones(size(xt));
            end

            uu_Big = max(max(abs(uu)));
            axes(handles.surfacePlots)
            surf(xt,tx,uu);
            xlabel('position'); ylabel('time'); zlabel('temperature');
            set(handles.upperPlot,'Value',2);
        end
    end

elseif module =='1DWave'
    [u,A_coefficients,B_coefficients,A_basis,B_basis,playString] = vibratingString(numModes, k,L,  IC_u, IC_dudt, alpha, beta);
    x_coords = linspace(0,L); t_coords = linspace(0,timeSpan,301);
    [xt,tx] = meshgrid(x_coords, t_coords);

    if diff(u,x) ~= 0
        uu = subs(u,{x,t}, {xt,tx});
    else
        uu = double(u) * ones(size(xt));
    end


    uu_Big = max(max(abs(uu)));
    uPlay = subs(playString,t,t_coords);
    axes(handles.surfacePlots)
    surf(xt,tx,uu);
    xlabel('position'); ylabel('time'); zlabel('displacement');
    set(handles.upperPlot,'Value',2);

elseif module == 'DmpStr'
    [u, A_coefficients, B_coefficients, A_basis, B_basis] = dampedString(numModes, k,L,IC_u, IC_dudt, alpha, beta, b);
    x_coords = linspace(0,L); t_coords = linspace(0,timeSpan,501);
    [xt,tx] = meshgrid(x_coords, t_coords);
    if diff(u,x) ~= 0
        uu = subs(u,{x,t}, {xt,tx});
    else
        uu = double(u) * ones(size(xt));
    end
    uu_Big = max(max(abs(uu)));
    axes(handles.surfacePlots)
    surf(xt,tx,uu);
    xlabel('position'); ylabel('time'); zlabel('displacement');
    set(handles.upperPlot,'Value',2);

elseif module == 'VbRect'
    set(handles.upperPlot,'Value',2);
    set(handles.pushbutton4,'Visible','Off');
    set(handles.edit9,'Visible','On');

    set(handles.edit9,'String','0');
    set(handles.slider2,'Visible','On');
    set(handles.slider2,'Value',0);
    set(handles.slider2,'Max',timeSpan);
    set(handles.slider2,'SliderStep',[1/100, 1/10]);

    [u, A_coefficients, B_coefficients, A_basis, B_basis] = vibratingRectangle(numModes, k, L, H, IC_u, IC_dudt,alpha,beta);
    wb = waitbar(0,'Evaluating . . .');
    x_coords = linspace(0,L,50); y_coords = linspace(0,H,50); t_coords = linspace(0,timeSpan,101);
    [xx,yy,tt] = ndgrid(x_coords, y_coords, t_coords);
    if and(diff(u,x) ~=0, diff(u,y) ~=0)
        uu = subs(u,{x,y,t},{xx,yy,tt});
    elseif diff(u,x)~=0
        uu = subs(u,{x,t},{xx,tt});
    elseif diff(u,y)~=0
        uu = subs(u,{y,t},{yy,tt});
    else
        uu = double(u) * ones(size(yy));
    end
    close(wb);
    uu_Big = max(abs(max(max(max(uu)))),abs(min(min(min(uu)))));
    axes(handles.surfacePlots)

    surf(uu(:,:,1)');
    set(handles.surfacePlots,'ZLim',[-uu_Big uu_Big]); set(handles.surfacePlots,'ZLimMode','Manual');
    xlabel('x position'); ylabel('y position'); zlabel('displacement');
elseif module == '2DHEAT'
    set(handles.upperPlot,'Value',2);
    set(handles.pushbutton4,'Visible','Off');

    set(handles.edit9,'Visible','On');

    set(handles.edit9,'String','0');
    set(handles.slider2,'Visible','On');
    set(handles.slider2,'Max',timeSpan);
    set(handles.slider2,'Value',0);
    set(handles.slider2,'SliderStep',[1/100, 1/10]);

    [u, A_coefficients, B_coefficients, A_basis, B_basis] = HeatRectangle(numModes, k, L, H, IC_u,alpha,beta);
    wb = waitbar(0,'Evaluating . . .');
    x_coords = linspace(0,L,50); y_coords = linspace(0,H,50); t_coords = linspace(0,timeSpan,101);
    [xx,yy,tt] = ndgrid(x_coords, y_coords, t_coords);

    if and(diff(u,x) ~=0, diff(u,y) ~=0)
        uu = subs(u,{x,y,t},{xx,yy,tt});
    elseif diff(u,x)~=0
        uu = subs(u,{x,t},{xx,tt});
    elseif diff(u,y)~=0
        uu = subs(u,{y,t},{yy,tt});
    else
        uu = double(u) * ones(size(yy));
    end

    close(wb);
    uu_Big = max(abs(max(max(max(uu)))),abs(min(min(min(uu)))));
    axes(handles.surfacePlots)

    surf(uu(:,:,1)');
    set(handles.surfacePlots,'ZLim',[-uu_Big uu_Big]); set(handles.surfacePlots,'ZLimMode','Manual');
    xlabel('x position'); ylabel('y position'); zlabel('temperature');

elseif module == '1DChar'
    cChar = get(handles.edit18,'String');

    [T,X_U,tt] = methodOfChars(cChar, source,numModes,timeSpan, IC_u, L, H);
    save TY T X_U
    uu_Big = max(abs(max(max(max(X_U(:,:,2))))));

    axes(handles.surfacePlots)
    save xu X_U
    mesh(X_U(1:tt,:,1),T,X_U(1:tt,:,2));
    xlabel('x position'); ylabel('time '); zlabel('u');
    set(handles.popupmenu1,'Value',6);
    set(handles.upperPlot,'Value',2);
end

if audioON
    uuPlay = norm(uu(1,:))*uPlay;
    save uPlay uuPlay
    sound(norm(uu(:,1)) * uPlay,500/2);
end

axes(handles.surfacePlots)
rotate3d;
set(handles.surfacePlots,'XMinorTick','on')
grid on
close(evalBar);



contents = get(handles.popupmenu1,'String');
selection = char(contents(get(handles.popupmenu1,'Value')));
set(handles.slider1,'Value',0);
set(handles.varvalue,'String','0');
lowerPlot(selection, handles);
axes(handles.surfacePlots)
rotate3d on;

function lowerPlot (whichPlot, handles)
hold off;
syms x t y;
global uu u x_coords t_coords y_coords xt tx xx yy tt ICF numModes timeSpan L H moduleType module modeSurf_subbed mS_Big;
global A_coefficients B_coefficients A_basis B_basis uu_Big tt;
global T X_U;
if strcmp(whichPlot,'fixed x') == 1
    axes(handles.linePlots)

    if moduleType == '1D'
        xValue = str2double(get(handles.xField, 'String'))
        if ~xValue
            set(handles.xField,'String',0);
            xValue = 0;
        end
        set(handles.xSlider,'Value',xValue);
        set(handles.xSlider,'SliderStep',[1/50 1/10]);

        set(handles.tSlider,'Visible','Off');
        set(handles.tField,'Visible','Off');
        set(handles.slider1,'Visible','Off');

        u_xSlice = subs(u,{x,t},{xValue * ones(1,length(t_coords)),t_coords});
        plot(t_coords, u_xSlice);
        xlabel('time'); ylabel('temperature');
        set(handles.linePlots,'YLim',[-uu_Big uu_Big]);
        set(handles.linePlots,'XMinorTick','on');
        grid on
    elseif moduleType == '2D'

        xValue = str2double(get(handles.xField,'String'))
        if xValue == NaN
            set(handles.xField,'String',0);
            xValue = 0;
        end

        set(handles.xSlider,'Value',xValue);

        xIndex = floor(49*xValue/L)+1;
        u_xSlice = uu(xIndex,:,:);

        axes(handles.linePlots);
        surf(squeeze(yy(xIndex,:,:)),squeeze(tt(xIndex,:,:)),squeeze(u_xSlice));
        set(handles.linePlots,'ZLim',[-uu_Big,uu_Big]); set(handles.linePlots,'ZLimMode','Manual');

        rotate3d;
    end

elseif strcmp(whichPlot,'fixed t') == 1
    axes(handles.linePlots)
    tValue = str2double(get(handles.tField, 'String'));


    if tValue == NaN
        set(handles.tField,'String',0);
        tValue = 0;
    end


    set(handles.tSlider,'Value',tValue);

    if moduleType == '1D'

        if module == '1DChar'
            'here'
            u_tSlice = X_U(1+floor(1000*tValue/timeSpan),:,2);
            x_coords = X_U(1+floor(1000*tValue/timeSpan),:,1);
        else
            u_tSlice = subs(u,{x,t},{x_coords,tValue*ones(1,length(x_coords))});
        end


        plot(x_coords, u_tSlice)
        set(handles.linePlots,'XMinorTick','on')
        xlabel('position'); ylabel('temperature');
        set(handles.linePlots,'YLim',[-1.2*uu_Big 1.2*uu_Big]);
        grid on
    elseif moduleType == '2D'
        t_index = floor(100*tValue/timeSpan)+1;
        u_tSlice = squeeze(uu(:,:,t_index));
        axes(handles.linePlots);
        surf(squeeze(xx(:,:,t_index)),squeeze(yy(:,:,t_index)),u_tSlice);
        set(handles.linePlots,'ZLim',[-uu_Big,uu_Big]); set(handles.linePlots,'ZLimMode','Manual');

    end

elseif strcmp(whichPlot,'Fourier Spectrum Magnitude') == 1
    set(handles.varvalue,'Visible','Off');
    set(handles.slider1,'Visible','Off');
    set(handles.edit13,'Visible','Off');
    axes(handles.linePlots);
    hold off;
    if strcmp(moduleType,'2D')
        bar3(.125:+1:(.125 + numModes), abs(A_coefficients),0.25,'b');
        pause
        hold on; bar3(.125:+1:(.125 + numModes), abs(B_coefficients),0.25,'r')
        set(handles.linePlots, 'ZLim',[0 1.2*max(max(max(abs(A_coefficients)),max(max(abs(B_coefficients)))))]);
        set(handles.linePlots,'XLim',[-.125 numModes+.25]);
        set(handles.linePlots,'YLim',[-.125 numModes+.25]);

        set(handles.linePlots,'XTick',[0:+1:numModes]);
        set(handles.linePlots,'YTick',[0:+1:numModes]);
        rotate3d;

    else
        bar(-.125:+1:(-.125 + numModes), abs(A_coefficients),0.25,'b');
        hold on; bar(.125:+1:(.125 + numModes), abs(B_coefficients),0.25,'r')
        max(max(abs(A_coefficients)),max(abs(B_coefficients)))
        set(handles.linePlots,'YLim',[0 1.2*max(max(abs(A_coefficients)),max(abs(B_coefficients)))]);
        set(handles.linePlots,'XTick',[0:+1:numModes])
        rotate3d;
    end
    grid on;
    hold off;
elseif strcmp(whichPlot,'Individual Fourier Modes') == 1
    set(handles.xField,'Visible','Off');
    set(handles.xSlider,'Visible','Off');

    if strcmp(moduleType,'2D')
        set(handles.slider1,'Visible','On');
        set(handles.slider1,'Max',timeSpan);
        set(handles.slider1,'SliderStep',[1/50 1/50]);
        set(handles.slider1,'Value',0);
        set(handles.varvalue,'String',0);

        nm= (get(handles.edit13,'String'));
        f = find(nm==',');
        n = str2double(nm(1:f-1)); m = str2double(nm(f+1:end));

        modeSurf = A_coefficients(n+1,m+1) * A_basis(n+1,m+1) + B_coefficients(n+1,m+1) * B_basis(n+1,m+1);
        wb = waitbar(1,'Evaluating . . .');


        if and(diff(modeSurf,x)~=0, diff(modeSurf,y) ~=0)
            modeSurf_subbed = subs(modeSurf,{x,y,t},{xx,yy,tt});
        elseif diff(modeSurf,x)~=0
            modeSurf_subbed = subs(modeSurf,{x,t},{xx,tt});
        elseif diff(modeSurf,y) ~= 0
            modeSurf_subbed = subs(modeSurf, {y,t},{yy,tt});
        else
            modeSurf_subbed = 0;
        end

        mS_Big = max(max(max(abs(modeSurf_subbed))));
        close(wb);
        axes(handles.linePlots);
        if and(max(size(modeSurf_subbed)) ~= 1, ~isempty(modeSurf_subbed))

            
            surf(modeSurf_subbed(:,:,1)');
            set(handles.linePlots,'ZLim',[-mS_Big,mS_Big]); set(handles.linePlots,'ZLimMode','Manual');
            xlabel('x position'); ylabel('y_position'); zlabel('displacement');
        else
            if or(modeSurf_subbed == 0, modeSurf_subbed == [])
                surf(zeros(100,100));
            else
                surf(double(modeSurf) * ones(100,100));
            end
        end

    elseif strcmp(moduleType,'1D')
        set(handles.varvalue,'Visible','On');

        nValue = str2double(get(handles.varvalue, 'String'));
        if ~nValue
            set(handles.varvalue,'String', 1);
        end
        n = str2double(get(handles.varvalue,'String'))
        set(handles.slider1,'Min',0);
        set(handles.slider1,'Max',numModes);
        set(handles.slider1, 'SliderStep',[1/numModes 1/numModes])

        modePlot = A_coefficients(n+1) * A_basis(n+1) + B_coefficients(n+1) * B_basis(n+1);
        wb = waitbar(0,'Evaluating . . .');
        modePlot_subbed = subs(modePlot, {x,t}, {xt,tx});
        close(wb);
        if max(size(modePlot_subbed)) ~= 1
            axes(handles.linePlots);
            surf(xt,tx,modePlot_subbed);
            rotate3d;
        else
            axes(handles.linePlots);
            surf(xt,tx,zeros(max(size(xt(:,1))),max(size(tx(1,:)))));
            rotate3d;
        end
    end

elseif strcmp(whichPlot, 'Characteristic Curves')

    axes(handles.linePlots);
    plot(X_U(1:tt,:,1),T);
    grid on;
    rotate3d off;
end


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in lower plot menu.
function popupmenu1_Callback(hObject, eventdata, handles)
global selection moduleType timeSpan L numModes;
contents = get(hObject,'String');
selection = char(contents(get(hObject,'Value')));

if strcmp(selection, 'fixed x') == 1
    hold off;
    set(handles.varvalue,'Visible','Off');
    set(handles.slider1,'Visible','Off');
    set(handles.edit13,'Visible','Off');

    set(handles.xSlider,'Visible','On');
    set(handles.xField,'Visible','On');

    set(handles.tSlider,'Visible','Off');
    set(handles.tField,'Visible','Off');

    set(handles.vtype,'Visible','On');

    set(handles.vtype, 'String','x');
    set(handles.xSlider,'Max',L);
    set(handles.xSlider,'SliderStep',[1/50 1/10]);
    lowerPlot ('fixed x',handles);
elseif strcmp(selection, 'fixed t') == 1
    set(handles.xSlider,'Visible','Off');
    set(handles.xField,'Visible','Off');

    set(handles.slider1,'Visible','Off');
    set(handles.varvalue,'Visible','Off');
    set(handles.edit13,'Visible','Off');

    set(handles.tField,'Visible','On');
    set(handles.tField,'String','0');
    set(handles.tSlider,'Visible','On');
    set(handles.tSlider,'Max',timeSpan);
    set(handles.tSlider,'SliderStep',[1/100 1/10]);

    set(handles.vtype,'Visible','On');
    set(handles.vtype, 'String','t');
    lowerPlot ('fixed t', handles);
elseif strcmp(selection, 'Fourier Spectrum Magnitude') == 1
    set(handles.varvalue,'Visible','Off');
    set(handles.slider1,'Visible','Off');

    set(handles.tSlider,'Visible','Off');
    set(handles.tField,'Visible','Off');

    set(handles.xField,'Visible','Off');
    set(handles.xSlider,'Visible','Off');

    set(handles.vtype,'Visible','Off');
    lowerPlot('Fourier Spectrum Magnitude',handles);
elseif strcmp(selection, 'Individual Fourier Modes') == 1
    set(handles.vtype,'Visible','Off');
    set(handles.varvalue,'Visible','On');

    set(handles.xSlider,'Visible','Off');
    set(handles.tSlider,'Visible','Off');
    set(handles.xField,'Visible','Off');
    set(handles.tField,'Visible','Off');

    set(handles.edit13,'Visible','On');

    set(handles.slider1,'Visible','On');
    set(handles.slider1,'Min',0);
    set(handles.slider1,'Max',numModes);
    set(handles.slider1,'SliderStep',[1/numModes 1/numModes]);


    if strcmp(moduleType,'2D')
        set(handles.slider1,'Max',timeSpan);
        set(handles.edit13,'String','1,1');
        set(handles.varvalue,'String','1');
    elseif strcmp(moduleType,'1D')
        set(handles.edit13,'Visible','Off');
        set(handles.slider1,'Visible','Off');
    end


    lowerPlot('Individual Fourier Modes',handles);

end



% --- Executes during object creation, after setting all properties.
function varvalue_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%% The mode numbers
function varvalue_Callback(hObject, eventdata, handles)
global selection;
set(handles.slider1,'Value',str2double(get(hObject,'String')));
lowerPlot(selection, handles);


function IC_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function IC_Callback(hObject, eventdata, handles)

function alpha_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function alpha_edit_Callback(hObject, eventdata, handles)
set(handles.popupmenu2,'Value',1);


function beta_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function beta_edit_Callback(hObject, eventdata, handles)
set(handles.popupmenu2,'Value',1);

function slider1_CreateFcn(hObject, eventdata, handles)
usewhitebg = 0;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
global selection L numModes timeSpan H modeSurf_subbed mS_Big moduleType;

sliderValue = get(hObject, 'Value');
if (strcmp(moduleType,'2D'))
    axes(handles.linePlots);
    surf(modeSurf_subbed(:,:,floor(100*sliderValue/timeSpan)+1));
    set(handles.linePlots,'ZLim',[-mS_Big,mS_Big]); set(handles.linePlots,'ZLimMode','Manual');
    xlabel('x position'); ylabel('y_position'); zlabel('displacement');
else
    lowerPlot(selection, handles);
end
try
    set(handles.varvalue, 'String', sliderValue);
end
% --- Executes during creation of Boundary Condition Menu
function popupmenu2_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in Boundary Condition Menu.
function popupmenu2_Callback(hObject, eventdata, handles)
BCSelection = get(hObject,'String');
BCSelectionValue = BCSelection(get(hObject,'Value'));
if strcmp(BCSelectionValue, 'Dirichlet')

    set(handles.alpha_edit,'String',1); set(handles.beta_edit,'String',0);
elseif strcmp(BCSelectionValue, 'Neumann')
    set(handles.alpha_edit,'String',0); set(handles.beta_edit,'String',1);
end






% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in Module Menu.
function popupmenu3_Callback(hObject, eventdata, handles)


global module audioON moduleType;
moduleSelectorContents = get(hObject,'String');
moduleSelection = moduleSelectorContents(get(hObject,'Value'));
global X_U uu pivotMovie;
clear global pivotMovie


switch char(moduleSelection)
    case '1-D Heat Equation'
        module = '1DHEAT';
        moduleType = '1D';
        audioON = 0;
    case '2-D Heat Equation'
        module = '2DHEAT';
        moduleType = '2D';
    case 'Vibrating String'
        module = '1DWave';
        moduleType = '1D';
    case 'Vibrating Rectangle'
        module = 'VbRect';
        moduleType = '2D';
        set(handles.pushbutton4,'Visible','Off');
        set(handles.edit9,'Visible','Off');

        set(handles.upperPlot,'Value',1);
    case 'Damped String'
        module = 'DmpStr'
        moduleType = '1D';
        set(handles.pushbutton4,'Visible','Off');
        set(handles.edit9,'Visible','Off');
        set(handles.upperPlot,'Value',1);
    case 'Method of Characteristics'
        module = '1DChar';
        moduleType = '1D';
        set(handles.pushbutton4,'Visible','Off');
        set(handles.edit9,'Visible','Off');

        set(handles.upperPlot,'Value',1);
end

eqWindow = eqInfo65;
movegui(eqWindow,'northwest');
fieldSetter(module,handles)
plotReset(handles);

% --- Executes on button press in audioButton.
function audioButton_Callback(hObject, eventdata, handles)
global audioON
if get(hObject,'Value') == 1
    audioON = 1
elseif get(hObject,'Value')==1
    audioON = 0
end

% --- Executes on button press in IC_plotter.
function IC_plotter_Callback(hObject, eventdata, handles)
global userIC_u userIC_dudt moduleType;
global leftEnd rightEnd
leftEnd = str2double(get(handles.L_input,'String'));
try
    rightEnd = str2double(get(handles.Height,'String'));
catch
    rightEnd = 0;
end

if moduleType  == '1D'
    [userIC_u, userIC_dudt] = plotIC;
    if length(userIC_u) ~= 0
        set(handles.IC,'String','USER-DEFINED');
    end

    if length(userIC_dudt) ~= 0
        set(handles.IC_derivative,'String','USER-DEFINED');
    end
else
    warndlg('The IC Plotter is currently only available for the 1 spatial dimension modules');
end
% --- Executes during creation of Upper Plot menu.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in Upper Plot Menu
function popupmenu4_Callback(hObject, eventdata, handles)
global module upperSelection uu xx yy moduleType xt tx X_U T timeSpan;
global pivotMovie;
contents = get(hObject,'String');
upperSelection = char(contents(get(hObject,'Value')))

if (strcmp(upperSelection, 'Movie'))
    set(handles.slider2,'Visible','Off');
    set(handles.edit9, 'Visible','Off');
    pivotMovie = movieMaker(handles);

    set(handles.pushbutton4,'Visible','On');

elseif strcmp(upperSelection, 'Surface Plot') == 1
    clear pivotMovie;
    set(handles.pushbutton4,'Visible','Off');
    'h'
    if module == '1DChar'
        axes(handles.surfacePlots);
        surf(X_U(1:300,:,1),T,X_U(1:300,:,2));
        rotate3d;
        size(X_U)
    elseif strcmp(moduleType,'2D')
        
        set(handles.edit9,'Visible','On');
        set(handles.slider2,'Visible','On');
        set(handles.slider2,'Max',timeSpan);
        set(handles.pushbutton4, 'Visible','Off');
        axes(handles.surfacePlots);
        surf(xx(:,:,1),yy(:,:,1),uu(:,:,1));
        rotate3d;
        set(handles.edit9,'String','0');
        set(handles.slider2,'Value',0);
    else
        set(handles.edit9,'Visible','Off');
        set(handles.slider2,'Visible','Off');
        set(handles.pushbutton4,'Visible','Off');
        axes(handles.surfacePlots);
        surf(xt,tx,uu);
        rotate3d;
    end
elseif strcmp(upperSelection, 'Mesh Plot') == 1
    clear pivotMovie
    if module == '1DChar'
       size(X_U)
        axes(handles.surfacePlots);
        mesh(X_U(1:300,:,1),T,X_U(1:300,:,2));
        rotate3d;
    elseif strcmp(moduleType,'2D')
        set(handles.edit9,'Visible','On');
        set(handles.slider2,'Visible','On');
        set(handles.pushbutton4, 'Visible','Off');
        axes(handles.surfacePlots);
        mesh(xx(:,:,1),yy(:,:,1),uu(:,:,1));
        set(handles.edit9,'String','0');
        set(handles.slider2,'Value',0);
        rotate3d;
    else
        set(handles.edit9,'Visible','Off');
        set(handles.slider2,'Visible','Off');
        set(handles.pushbutton4,'Visible','Off');
        axes(handles.surfacePlots);
        mesh(xt,tx,uu);
        rotate3d;
    end
elseif strcmp(upperSelection,'Top View')
    clear pivotMovie
    axes(handles.surfacePlots)
    if module == '1DChar'

        pcolor(real(X_U(1:300,:,1)),T,real(X_U(1:300,:,2))), shading interp
        grid on;
    elseif moduleType == '1D'
        pcolor(xt,tx,uu), shading interp;

    end
end
% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%% edit9 = upper time value
function edit9_Callback(hObject, eventdata, handles)
me = get(handles.edit9,'String')
set(handles.slider2,'Value',str2double(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
global timeSpan

usewhitebg = 1;
set(hObject,'Min',0);

if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on upper Slider Movement
function slider2_Callback(hObject, eventdata, handles)
global timeSpan uu uu_Big xx yy upperSelection;
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliderValue = get(hObject, 'Value');
set(handles.edit9, 'String', sliderValue);
time_index = floor(sliderValue/timeSpan * 100) + 1;

axes(handles.surfacePlots)

upperSelection
if strcmp(upperSelection, 'Surface Plot')
    surf(xx(:,:,1),yy(:,:,1),uu(:,:,time_index));
    rotate3d;
elseif strcmp(upperSelection, 'Mesh Plot')
    mesh(xx(:,:,1),yy(:,:,1),uu(:,:,time_index));
    rotate3d;
end

set(handles.surfacePlots,'ZLim',[-uu_Big uu_Big]); set(handles.surfacePlots,'ZLimMode','Manual');
xlabel('x position'); ylabel('y_position'); zlabel('displacement');

function thisMovie = movieMaker(handles)
global module uu uu_Big xx yy x_coords t_coords L H moduleType T X_U;
axes(handles.surfacePlots);

if ~strcmp(module,'1DChar')

    if strcmp(moduleType,'2D')
        for i = 1:max(size(uu(1,1,:)))
            %wb = waitbar(i/max(size(uu(1,1,:))),'Making Movie . . .');

            surf(xx(:,:,1),yy(:,:,1),uu(:,:,i));
            set(gca,'ZLim',[-uu_Big uu_Big]);
            thisMovie(i) = getframe;
        end
    else
        for i = 1:size(uu(:,1))
            %wb = waitbar(i/max(size(uu(1,1,:))),'Making Movie . . .');
            plot(x_coords,uu(i,:));
            set(gca,'YLim',[-uu_Big uu_Big]);
            grid on;
            thisMovie(i) = getframe;
        end

    end
    %     close(wb);
else
    for i = 1:+5:500
        plot(X_U(i,:,1),X_U(i,:,2));
        set(gca,'YLim',[-uu_Big uu_Big]);
        grid on;
        thisMovie(i) = getframe;
    end
    %close(wb);
end

% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('pivotHelp.html');

% --------------------------------------------------------------------
function HelpMenu_Callback(hObject, eventdata, handles)
% hObject    handle to HelpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Height_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Height_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function IC_derivative_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function IC_derivative_Callback(hObject, eventdata, handles)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
global pivotMovie
axes(handles.surfacePlots);
movie(pivotMovie,0);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aboutPivot;


% --------------------------------------------------------------------
function properties_Callback(hObject, eventdata, handles)
% hObject    handle to properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%% Mode number editable field
function edit13_Callback(hObject, eventdata, handles)
lowerPlot('Individual Fourier Modes',handles);


% --------------------------------------------------------------------
function autoScale_Callback(hObject, eventdata, handles)

% hObject    handle to lowerZoomOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.linePlots);
set(handles.linePlots,'ZLimMode','Auto');
xlabel('x position'); ylabel('y_position'); zlabel('displacement');

% --------------------------------------------------------------------
function resetScale_Callback(hObject, eventdata, handles)
global uu_Big;

% hObject    handle to lowerZoomOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.linePlots);
set(handles.linePlots,'ZLim',[-uu_Big,uu_Big]); set(handles.linePlots,'ZLimMode','Manual');
xlabel('x position'); ylabel('y_position'); zlabel('displacement');


% --------------------------------------------------------------------
function linePlotMenu_Callback(hObject, eventdata, handles)
% hObject    handle to linePlotMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function lpMenu_Callback(hObject, eventdata, handles)
% hObject    handle to lpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function b_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to b_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function b_input_Callback(hObject, eventdata, handles)
% hObject    handle to b_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of b_input as text
%        str2double(get(hObject,'String')) returns contents of b_input as a double


% --- Executes during object creation, after setting all properties.
function sourceFxn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sourceFxn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function sourceFxn_Callback(hObject, eventdata, handles)
% hObject    handle to sourceFxn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sourceFxn as text
%        str2double(get(hObject,'String')) returns contents of sourceFxn as a double


% --- Executes during object creation, after setting all properties.
function timeSpan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeSpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function timeSpan_Callback(hObject, eventdata, handles)
% hObject    handle to timeSpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeSpan as text
%        str2double(get(hObject,'String')) returns contents of timeSpan as a double


% --- Executes during object creation, after setting all properties.
function tSpan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tSpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function tSpan_Callback(hObject, eventdata, handles)
% hObject    handle to tSpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tSpan as text
%        str2double(get(hObject,'String')) returns contents of tSpan as a double


% --------------------------------------------------------------------
function zmON_Callback(hObject, eventdata, handles)
% hObject    handle to zmON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

get(hObject,'Checked')

if strcmp(get(hObject,'Checked'),'off')
    'here?'
    set(hObject,'Checked','On');
    axes(handles.surfacePlots);
    set(gca,'CameraViewAngle',get(gca,'CameraViewAngle')-5)
elseif strcmp(get(hObject,'Checked'),'on')
    set(hObject,'Checked','Off');
    axes(handles.surfacePlots);
    zoom off;
end

% --------------------------------------------------------------------
function zmOFF_Callback(hObject, eventdata, handles)
% hObject    handle to zmOFF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function fieldSetter(module,handles)
if module == '1DHEAT'
    set(handles.f2_input,'Visible','On');
    set(handles.text19, 'Visible','Off');
    set(handles.IC_derivative,'Visible','Off');
    set(handles.Height,'Visible','Off');
    set(handles.text_height,'Visible','Off');
    set(handles.b_input,'Visible','Off');
    set(handles.text16,'Visible','Off');
    set(handles.text3,'String','k');
    set(handles.edit18,'Visible','Off');
    set(handles.c_forChars,'Visible','Off');
    set(handles.alpha_edit,'Visible','On');
    set(handles.beta_edit,'Visible','On');
    set(handles.popupmenu2,'Visible','On');
elseif module == '1DWave'
    set(handles.f2_input,'Visible','On');
    set(handles.text19, 'Visible','On');
    set(handles.IC_derivative,'Visible','On');
    set(handles.Height,'Visible','Off');
    set(handles.text_height,'Visible','Off');
    set(handles.b_input,'Visible','Off');
    set(handles.text16,'Visible','Off');
    set(handles.text3,'String','c');
    set(handles.edit18,'Visible','Off');
    set(handles.c_forChars,'Visible','Off');
    set(handles.alpha_edit,'Visible','On');
    set(handles.beta_edit,'Visible','On');
    set(handles.popupmenu2,'Visible','On');
elseif module == '2DHEAT'

    set(handles.f2_input,'Visible','On');
    set(handles.text_height,'String','H');
    set(handles.text19, 'Visible','Off');
    set(handles.IC_derivative,'Visible','Off');
    set(handles.Height,'Visible','On');
    set(handles.text_height,'Visible','On');
    set(handles.b_input,'Visible','Off');
    set(handles.text16,'Visible','Off');
    set(handles.text3,'String','k');
    set(handles.edit18,'Visible','Off');
    set(handles.c_forChars,'Visible','Off');
    set(handles.alpha_edit,'Visible','On');
    set(handles.beta_edit,'Visible','On');
    set(handles.popupmenu2,'Visible','On');


elseif module == 'VbRect'

    set(handles.f2_input,'Visible','On');
    set(handles.text_height,'String','H');
    set(handles.text19, 'Visible','On');
    set(handles.IC_derivative,'Visible','On');
    set(handles.Height,'Visible','On');
    set(handles.text_height,'Visible','On');
    set(handles.b_input,'Visible','Off');
    set(handles.text16,'Visible','Off');
    set(handles.text3,'String','c');
    set(handles.edit18,'Visible','Off');
    set(handles.c_forChars,'Visible','Off');
    set(handles.alpha_edit,'Visible','On');
    set(handles.beta_edit,'Visible','On');
    set(handles.popupmenu2,'Visible','On');
elseif module == 'DmpStr'
    set(handles.f2_input,'Visible','On');
    set(handles.text19, 'Visible','On');
    set(handles.IC_derivative,'Visible','On');
    set(handles.Height,'Visible','Off');
    set(handles.text_height,'Visible','Off');
    set(handles.b_input,'Visible','On');
    set(handles.text16,'Visible','On');
    set(handles.text3,'String','c');
    set(handles.edit18,'Visible','Off');
    set(handles.c_forChars,'Visible','Off');
    set(handles.alpha_edit,'Visible','On');
    set(handles.beta_edit,'Visible','On');
    set(handles.popupmenu2,'Visible','On');

elseif module == '1DChar'
    set(handles.f2_input,'Visible','Off');
    set(handles.text19, 'Visible','On');
    set(handles.IC_derivative,'Visible','Off');
    set(handles.Height,'Visible','On');
    set(handles.text_height,'Visible','On');
    set(handles.text_height,'String','R');
    set(handles.b_input,'Visible','Off');
    set(handles.text16,'Visible','Off');
    set(handles.text3,'String','c');
    set(handles.edit18,'Visible','On');
    set(handles.c_forChars,'Visible','On');
    set(handles.alpha_edit,'Visible','Off');
    set(handles.beta_edit,'Visible','Off');
    set(handles.popupmenu2,'Visible','Off');

end





% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double

function plotReset(handles)
axes(handles.surfacePlots);
surf(zeros(100,100));
xlabel('position'),ylabel('time'),zlabel('temperature');
axes(handles.linePlots);
plot(zeros(1,10));
grid on;
set(handles.slider1,'Visible','Off');
set(handles.slider2,'Visible','Off');
set(handles.varvalue,'Visible','Off');
set(handles.edit13,'Visible','Off');


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global xLeft xRight yLeft yRight moduleType ;

warndlg('The Boundary Condition Editor to allow for non-Homogeneous BCs is still under construction');
% [xLeft, xRight, yLeft, yRight] = BCeditor;
% set(handles.alpha_edit,'String','Edited');
% set(handles.beta_edit,'String','Edited');

% --- Executes on button press in printButton.
function printButton_Callback(hObject, eventdata, handles)
% hObject    handle to printButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg;



function xField_Callback(hObject, eventdata, handles)
% hObject    handle to xField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xField as text
%        str2double(get(hObject,'String')) returns contents of xField as a double
global selection
set(handles.xSlider,'Value',str2double(get(hObject,'String')));
lowerPlot(selection, handles);

% --- Executes during object creation, after setting all properties.
function xField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function xSlider_Callback(hObject, eventdata, handles)
% hObject    handle to xSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sliderValue = get(hObject, 'Value');
set(handles.xField,'String',sliderValue);
lowerPlot('fixed x',handles);

% --- Executes during object creation, after setting all properties.
function xSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on slider movement.
function tSlider_Callback(hObject, eventdata, handles)
% hObject    handle to tSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderValue = get(hObject, 'Value');
set(handles.tField,'String',sliderValue);
lowerPlot('fixed t',handles);

% --- Executes during object creation, after setting all properties.
function tSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function tField_Callback(hObject, eventdata, handles)
% hObject    handle to tField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tField as text
%        str2double(get(hObject,'String')) returns contents of tField as a double


% --- Executes during object creation, after setting all properties.
function tField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function saveState_Callback(hObject, eventdata, handles)
% hObject    handle to saveState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


