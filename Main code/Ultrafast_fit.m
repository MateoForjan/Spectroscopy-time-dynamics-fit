function varargout = Ultrafast_fit(varargin)
% ULTRAFAST_FIT MATLAB code for Ultrafast_fit.fig


%--------------------------------------------------------------------------
% Code begins with defining functions for loading a file, plotting and
% fitting.
% Callbacks for all the buttons follow up.
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% DESCRIPTION:

% Tool used for fitting functions on time dynamics obtained by measurements.
% 
% Singular value decomposition (SVD) analysis is automatically done at start. 
% Singular values and vectors (spectra and time dynamics) can be explored 
% by choosing the number of SVD values. 
% SVD analysis is of great help when considering doing global analysis.
% 
% Function for fitting: convolution of exponential function 
% (decay or rising dynamics) 
% with a gaussian (representing instrument response function). Additionally, 
% a gaussian and first two derivatives are summed with mentioned convolution 
% for the description of coherent artifact near time zero.
% 
% Residuals of fitting and subtracted coherent artifact can be explored.
% 
% Details are given in MANUAL.
% 
% -Choosing the 'wavelength to fit' plots a time-cut at chosen wavelength.
% 
% -Changing 'Nr. of exponentials' and initial coefficients plots the initial 
% guess for a fit.
% 
% -'Fit Curve' fits the curve and in pop-up windows gives the fit results 
% and goodness of fit.
%--------------------------------------------------------------------------



% Last Modified by GUIDE v2.5 25-Apr-2023 12:54:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Ultrafast_fit_OpeningFcn, ...
                   'gui_OutputFcn',  @Ultrafast_fit_OutputFcn, ...
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


% --- Executes just before Ultrafast_fit is made visible.
function Ultrafast_fit_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Ultrafast_fit (see VARARGIN)

% Choose default command line output for Ultrafast_fit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Ultrafast_fit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Ultrafast_fit_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%--------------------------------------------------------------------------
%           FUNCTIONS FOR LOADING FILES, PLOTTING AND FITTING
%--------------------------------------------------------------------------



% --- LOADING A FILE (load file button)
function loadfile_Callback(~, ~, handles)

global Wavelength
global Time
global Transient
global fajl
global clrbar
global axesfont
global fontsize
global checker_fit

checker_fit = 0; %checker used later;

fajl = uigetfile('*.dat');
timefajl = uigetfile('*.dat');
lambdafajl = uigetfile('*.dat');
Transient=dlmread(fajl);
Time=dlmread(timefajl);
Wavelength=dlmread(lambdafajl);
Time(end) = [];
Transient(:,end) = [];
Transient(isnan(Transient)) = 0;
Transient(isinf(Transient)) = 0;


%SETTING INITIAL VALUES
current_folder = pwd;
set(handles.statictext1, 'String', current_folder);

current_folder = fajl;
set(handles.static2, 'String', fajl);

fontsize=10;
clrbar=10;
axesfont = 10;

plot2d(handles);
plot_cut(Wavelength(round(end/2)),handles);

set(handles.Wavelength_fit, 'String', num2str(Wavelength(round(end/2))));

% PLOT OF SVD VALUES
global U
global S
global V

axes(handles.axes3);
cla(handles.axes3);
view(3)
[U,S,V] = svd(Transient);
SVD = diag(S);
g1=plot(SVD,'-o');set(g1,'LineWidth',4);hold on

line(xlim(), [0,0], 'LineWidth', 2, 'Color', 'k');
axis([-inf 10 -inf inf]); %340 650
xlabel(['SVD number'],'FontSize',12);
ylabel(['SVD value'],'FontSize',12);
% yticks([-1 -0.5 0 0.5 1 1.5 2])
xt = get(gca, 'XTick');set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperPosition', [0 0 15 15]);set(gca, 'FontSize', 12);box on;grid on;


plot_svd_dynamics(1,handles);
plot_svd_spectrum(1,handles);


%FUNCTION FOR PLOTTING A MATRIX
function plot2d(handles)
global Wavelength
global Time
global Transient
global t1
global t2
global lambda_1
global lambda_2
global figure1
global C2
global C1
global clrbar
global axesfont
global fontsize

axes(handles.axes2);
cla(handles.axes2);
figure1 = pcolor(handles.axes2,Time,Wavelength,Transient);
view(2)
colormap jet;
set(gca,'TickLength',[0.03 0.03]);
C1 =str2double(get(handles.C1,'String'));
C2 = str2double(get(handles.C2,'String'));
caxis([C1 C2]);
t1 = str2double(get(handles.Time_1,'String'));
t2 = str2double(get(handles.Time_2,'String'));
lambda_1 = str2double(get(handles.Wavelength1,'String'));
lambda_2 = str2double(get(handles.Wavelength2,'String'));
axis([t1 t2 lambda_1 lambda_2]);
set(gca,'XColor', 'k');
set(gca,'Layer','top');
set(gca,'FontSize',fontsize);
set(gca,'GridLineStyle','none');
colorbar('FontSize',clrbar);
shading interp;
xlabel('Delay [fs]','FontSize',axesfont,'Color','k');
ylabel('Wavelength (nm)','FontSize',axesfont,'Color','k');
addToolbarExplorationButtons(gcf)
hold on;
box on;


%FUNCTION FOR PLOTTING SVD DYNAMICS
function plot_svd_dynamics(r,handles)

global Time
global fontsize
global axesfont
global g1
global V
global t1
global t2


axes(handles.axes4);
% cla(handles.axes4);
x = 2;
y=2; %HOW MANY PIXELS TO SMOOTH
g1=plot(handles.axes4,Time,smooth(V(:,r),y),'DisplayName',num2str(r));set([g1],'LineWidth',x);
grid on;
t1 = str2double(get(handles.Time_1,'String'));
t2 = str2double(get(handles.Time_2,'String'));
lambda_1 = str2double(get(handles.Wavelength1,'String'));
lambda_2 = str2double(get(handles.Wavelength2,'String'));
axis([t1 t2 -inf inf]);
addToolbarExplorationButtons(gcf)
lgd = legend();
xlabel(['Delay [fs]'],'FontSize',axesfont);
ylabel('a.u.','FontSize',axesfont);
xt = get(gca, 'XTick');set(gca, 'FontSize', fontsize);box on;grid on;

%FUNCTION FOR PLOTTING SVD SPECTRA
function plot_svd_spectrum(r,handles)

global Wavelength
global fontsize
global axesfont
global h1
global U
global lambda_1
global lambda_2

x = 2;
y=2;
axes(handles.axes5);
h1=plot(handles.axes5,Wavelength,smooth(U(:,r),y),'DisplayName',[num2str(r)]);set([h1],'LineWidth',x);
grid on;
lambda_1 = str2double(get(handles.Wavelength1,'String'));
lambda_2 = str2double(get(handles.Wavelength2,'String'));
axis([lambda_1 lambda_2 -inf inf ]);
lgd = legend();
xlabel(['Wavelength [nm]'],'FontSize',axesfont);
ylabel('a.u.','FontSize',axesfont);
xt = get(gca, 'XTick');set(gca, 'FontSize', fontsize);box on;grid on;

%FUNCTION FOR PLOTTING THE DYNAMICS FOR FITTING
function plot_cut(valna1,handles)

global Wavelength
global Time
global Transient
global t1
global t2
global fontsize
global axesfont
global h2
global dynamics


t1 =str2double(get(handles.Time_1,'String'));
t2 =str2double(get(handles.Time_2,'String'));
I = find(abs(Wavelength-valna1)<((Wavelength(end)-Wavelength(end-1))/2));a=round(size(I)/2);a = a(1);I = I(a);Wavelength_I = Wavelength(I);
axes(handles.axes6);
% cla(handles.axes2);
delete(h2);
h2=plot(handles.axes6,Time,Transient(I,:),'color',[00 0.4 1],'DisplayName',[num2str(round(valna1)) ' nm']);set([h2],'LineWidth',2)
dynamics = Transient(I,:);
grid on;
axis([t1 t2 -inf inf]);
addToolbarExplorationButtons(gcf)
lgd = legend();

xlabel(['Delay [fs]'],'FontSize',axesfont);
ylabel('\DeltaA','FontSize',axesfont);
xt = get(gca, 'XTick');set(gca, 'FontSize', fontsize);box on;grid on;

%FUNCTION FOR PLOTTING THE INITIAL GUESS
function plot_ansatz(~,handles)

global Time
global A1
global A2
global A3
global s
global c
global B1
global B2
global B3
global B4
global tau1
global tau2
global tau3
global tau4
global h3
global t1
global t2
global h4

%POP = NUMBER OF EXPONENTIALS
pop = get(handles.popupmenu3,'Value');

if pop == 1
    delete(h3);
    axes(handles.axes6);
    delete(h4);
    h3 = plot(handles.axes6, Time, A1*exp(-(Time-c).^2 /(2*s^2)) + A2*(c-Time)/s^3.*exp(-(Time-c).^2/(2*s^2)) + A3*(Time.^2-2*c*Time-s^2 + c^2)/s^5.*exp(-(Time-c).^2 /(2*s^2)) +B1*exp(-Time/tau1).*exp( (c+s^2/(2*tau1))/tau1).*(1+erf( (Time-c-(s^2/tau1))/(1.41421356237*s))) , 'r');
    axis([t1 t2 -0.001 inf]);
    set([h3],'LineWidth',2);  
    grid on;
elseif pop == 2
    delete(h3);
    axes(handles.axes6);
    delete(h4);
    h3 = plot(handles.axes6, Time, A1*exp(-(Time-c).^2 /(2*s^2)) + A2*(c-Time)/s^3.*exp(-(Time-c).^2/(2*s^2)) + A3*(Time.^2-2*c*Time-s^2 + c^2)/s^5.*exp(-(Time-c).^2 /(2*s^2)) +B1*exp(-Time/tau1).*exp( (c+s^2/(2*tau1))/tau1).*(1+erf( (Time-c-(s^2/tau1))/(1.41421356237*s)))+B2*exp(-Time/tau2).*exp( (c+s^2/(2*tau2))/tau2).*(1+erf( (Time-c-(s^2/tau2))/(1.41421356237*s))) , 'r');
    axis([t1 t2 -0.001 inf]);
    set([h3],'LineWidth',2);
    grid on;
elseif pop == 3
    delete(h3);
    axes(handles.axes6);
    delete(h4);
    h3 = plot(handles.axes6, Time, A1*exp(-(Time-c).^2 /(2*s^2)) + A2*(c-Time)/s^3.*exp(-(Time-c).^2/(2*s^2)) + A3*(Time.^2-2*c*Time-s^2 + c^2)/s^5.*exp(-(Time-c).^2 /(2*s^2)) +B1*exp(-Time/tau1).*exp( (c+s^2/(2*tau1))/tau1).*(1+erf( (Time-c-(s^2/tau1))/(1.41421356237*s)))+B2*exp(-Time/tau2).*exp( (c+s^2/(2*tau2))/tau2).*(1+erf( (Time-c-(s^2/tau2))/(1.41421356237*s)))+B3*exp(-Time/tau3).*exp( (c+s^2/(2*tau3))/tau3).*(1+erf( (Time-c-(s^2/tau3))/(1.41421356237*s))) , 'r');
    axis([t1 t2 -0.001 inf]);
    set([h3],'LineWidth',2);
    grid on;
elseif pop == 4
    delete(h3);
    axes(handles.axes6);
    delete(h4);
    h3 = plot(handles.axes6, Time, A1*exp(-(Time-c).^2 /(2*s^2)) + A2*(c-Time)/s^3.*exp(-(Time-c).^2/(2*s^2)) + A3*(Time.^2-2*c*Time-s^2 + c^2)/s^5.*exp(-(Time-c).^2 /(2*s^2)) +B1*exp(-Time/tau1).*exp( (c+s^2/(2*tau1))/tau1).*(1+erf( (Time-c-(s^2/tau1))/(1.41421356237*s)))+B2*exp(-Time/tau2).*exp( (c+s^2/(2*tau2))/tau2).*(1+erf( (Time-c-(s^2/tau2))/(1.41421356237*s)))+B3*exp(-Time/tau3).*exp( (c+s^2/(2*tau3))/tau3).*(1+erf( (Time-c-(s^2/tau3))/(1.41421356237*s)))+B4*exp(-Time/tau4).*exp( (c+s^2/(2*tau4))/tau4).*(1+erf( (Time-c-(s^2/tau4))/(1.41421356237*s))) , 'r');
    axis([t1 t2 -0.001 inf]);
    set([h3],'LineWidth',2);  
    grid on;
end


%FUNCTION FOR PLOTTING A FIT
function plot_fit(~,handles)
global Time
global h3
global t1
global t2
global h4
global coef_A1
global coef_A2
global coef_A3
global coef_B1
global coef_c
global coef_s
global coef_tau1
global coef_B2
global coef_tau2


pop = get(handles.popupmenu3,'Value');

if pop == 1
    delete(h3);
    axes(handles.axes6);
    delete(h4);
    h3 = plot(handles.axes6, Time, coef_A1*exp(-(Time-coef_c).^2 /(2*coef_s^2)) + coef_A2*(coef_c-Time)/coef_s^3.*exp(-(Time-coef_c).^2/(2*coef_s^2)) + coef_A3*(Time.^2-2*coef_c*Time-coef_s^2 + coef_c^2)/coef_s^5.*exp(-(Time-coef_c).^2 /(2*coef_s^2)) +coef_B1*exp(-Time/coef_tau1).*exp( (coef_c+coef_s^2/(2*coef_tau1))/coef_tau1).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau1))/(1.41421356237*coef_s))),'DisplayName','Fit result','color',[1 0 0]  );
    axis([t1 t2 -0.001 inf]);
    set([h3],'LineWidth',2);  
    grid on;
elseif pop == 2
    global coef_B2
    global coef_tau2
    delete(h3);
    axes(handles.axes6);
    delete(h4);
    h3 = plot(handles.axes6, Time, coef_A1*exp(-(Time-coef_c).^2 /(2*coef_s^2)) + coef_A2*(coef_c-Time)/coef_s^3.*exp(-(Time-coef_c).^2/(2*coef_s^2)) + coef_A3*(Time.^2-2*coef_c*Time-coef_s^2 + coef_c^2)/coef_s^5.*exp(-(Time-coef_c).^2 /(2*coef_s^2)) +coef_B1*exp(-Time/coef_tau1).*exp( (coef_c+coef_s^2/(2*coef_tau1))/coef_tau1).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau1))/(1.41421356237*coef_s)))+coef_B2*exp(-Time/coef_tau2).*exp( (coef_c+coef_s^2/(2*coef_tau2))/coef_tau2).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau2))/(1.41421356237*coef_s))),'DisplayName','Fit result','color',[1 0 0]  );
    axis([t1 t2 -0.001 inf]);
    set([h3],'LineWidth',2);
    grid on;
    
elseif pop == 3
    global coef_B2
    global coef_tau2
    global coef_B3
    global coef_tau3
    delete(h3);
    axes(handles.axes6);
    delete(h4);
    h3 = plot(handles.axes6, Time, coef_A1*exp(-(Time-coef_c).^2 /(2*coef_s^2)) + coef_A2*(coef_c-Time)/coef_s^3.*exp(-(Time-coef_c).^2/(2*coef_s^2)) + coef_A3*(Time.^2-2*coef_c*Time-coef_s^2 + coef_c^2)/coef_s^5.*exp(-(Time-coef_c).^2 /(2*coef_s^2)) +coef_B1*exp(-Time/coef_tau1).*exp( (coef_c+coef_s^2/(2*coef_tau1))/coef_tau1).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau1))/(1.41421356237*coef_s)))+coef_B2*exp(-Time/coef_tau2).*exp( (coef_c+coef_s^2/(2*coef_tau2))/coef_tau2).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau2))/(1.41421356237*coef_s)))+coef_B3*exp(-Time/coef_tau3).*exp( (coef_c+coef_s^2/(2*coef_tau3))/coef_tau3).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau3))/(1.41421356237*coef_s))),'DisplayName','Fit result','color',[1 0 0]  );
    axis([t1 t2 -0.001 inf]);
    set([h3],'LineWidth',2);
    grid on;
elseif pop == 4
    global coef_B2
    global coef_tau2
    global coef_B3
    global coef_tau3
    global coef_B4
    global coef_tau4
    delete(h3);
    axes(handles.axes6);
    delete(h4);
    h3 = plot(handles.axes6, Time, coef_A1*exp(-(Time-coef_c).^2 /(2*coef_s^2)) + coef_A2*(coef_c-Time)/coef_s^3.*exp(-(Time-coef_c).^2/(2*coef_s^2)) + coef_A3*(Time.^2-2*coef_c*Time-coef_s^2 + coef_c^2)/coef_s^5.*exp(-(Time-coef_c).^2 /(2*coef_s^2)) +coef_B1*exp(-Time/coef_tau1).*exp( (coef_c+coef_s^2/(2*coef_tau1))/coef_tau1).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau1))/(1.41421356237*coef_s)))+coef_B2*exp(-Time/coef_tau2).*exp( (coef_c+coef_s^2/(2*coef_tau2))/coef_tau2).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau2))/(1.41421356237*coef_s)))+coef_B3*exp(-Time/coef_tau3).*exp( (coef_c+coef_s^2/(2*coef_tau3))/coef_tau3).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau3))/(1.41421356237*coef_s)))+coef_B4*exp(-Time/coef_tau4).*exp( (coef_c+coef_s^2/(2*coef_tau4))/coef_tau4).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau4))/(1.41421356237*coef_s))),'DisplayName','Fit result','color',[1 0 0]  );
    axis([t1 t2 -0.001 inf]);
    set([h3],'LineWidth',2);  
    grid on;
end

%FUNCTION FOR FITTING - POP = NUMBER OF EXPONENTIAL FUNCTIONS (chosen in
%pop up menu)
function fit_and_plot(pop, handles)

global dynamics
global Time
global A1
global A2
global A3
global B1
global B2
global B3
global B4
global c
global s
global tau1
global tau2
global tau3
global tau4
global t1 
global t2
global h4
global A1_min
global A2_min
global A3_min
global c_min
global s_min
global B1_min
global B2_min
global B3_min
global B4_min
global tau_1_min
global tau_2_min
global tau_3_min
global tau_4_min
global A1_max
global A2_max
global A3_max
global c_max
global s_max
global B1_max
global B2_max
global B3_max
global B4_max
global tau_1_max
global tau_2_max
global tau_3_max
global tau_4_max
global h3
global fit_line
global XPM
global coef_A1
global coef_A2
global coef_A3
global coef_B1
global coef_c
global coef_s
global coef_tau1
% global pop

global checker_fit

checker_fit = 1; %Checker used later

pop = get(handles.popupmenu3,'Value');

if pop == 1
    yfit1 = fittype('A1*exp(-(t-c)^2 /(2*s^2)) + A2*(c-t)/s^3* exp(-(t-c)^2 /(2*s^2))+ A3*(t^2-2*c*t-s^2 + c^2)/s^5*exp(-(t-c)^2 /(2*s^2)) + B1*exp(-t/tau1)*exp( (c+s^2/(2*tau1))/tau1)*(1+erf( (t-c-(s^2/tau1))/(1.41421356237*s)))' ,'ind','t'); 
    opts=fitoptions(yfit1);
    set(opts,'robust','on','Algorithm', 'Levenberg-Marquardt','tolfun',1e-16,'tolx',1e-16,'maxiter',10000,'MaxFunEvals',10000, ...
        'StartPoint',[A1 A2 A3 B1 c s tau1 ],'Lower', [A1_min A2_min A3_min B1_min c_min s_min tau_1_min],'Upper', [A1_max A2_max A3_max B1_max c_max s_max tau_1_max] ,'DiffMinChange', 1.0000e-08, 'DiffMaxChange', 10000);
    [fr , gf]= fit(Time,dynamics',yfit1,opts)

    coef_A1 = fr.A1;
    coef_A2 = fr.A2;
    coef_A3 = fr.A3;
    coef_B1 = fr.B1;
    coef_c = fr.c;
    coef_s = fr.s;
    coef_tau1 = fr.tau1;

    
    axes(handles.axes6);
    delete(h4);
    delete(h3);
    h4 = plot(handles.axes6, Time, coef_A1*exp(-(Time-coef_c).^2 /(2*coef_s^2)) + coef_A2*(coef_c-Time)/coef_s^3.*exp(-(Time-coef_c).^2/(2*coef_s^2)) + coef_A3*(Time.^2-2*coef_c*Time-coef_s^2 + coef_c^2)/coef_s^5.*exp(-(Time-coef_c).^2 /(2*coef_s^2)) +coef_B1*exp(-Time/coef_tau1).*exp( (coef_c+coef_s^2/(2*coef_tau1))/coef_tau1).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau1))/(1.41421356237*coef_s))),'DisplayName','Fit result','color',[1 0 0]  );
    axis([t1 t2 -inf inf]);
    set([h4],'LineWidth',2);  
    grid on;

elseif pop ==2
    yfit1 = fittype('A1*exp(-(t-c)^2 /(2*s^2)) + A2*(c-t)/s^3* exp(-(t-c)^2 /(2*s^2))+ A3*(t^2-2*c*t-s^2 + c^2)/s^5*exp(-(t-c)^2 /(2*s^2)) + B1*exp(-t/tau1)*exp( (c+s^2/(2*tau1))/tau1)*(1+erf( (t-c-(s^2/tau1))/(1.41421356237*s))) + B2*exp(-t/tau2)*exp( (c+s^2/(2*tau2))/tau2)*(1+erf( (t-c-(s^2/tau2))/(1.41421356237*s)))','ind','t'); 
    opts=fitoptions(yfit1);
    set(opts,'robust','off','Algorithm', 'Levenberg-Marquardt','tolfun',1e-18,'tolx',1e-18,'maxiter',10000,'MaxFunEvals',10000, ...
        'StartPoint',[A1 A2 A3 B1 B2 c s tau1 tau2 ],'Lower', [A1_min A2_min A3_min B1_min B2_min c_min s_min tau_1_min tau_2_min],'Upper', [A1_max A2_max A3_max B1_max B2_max c_max s_max tau_1_max tau_2_max], 'DiffMinChange', 1.0000e-08, 'DiffMaxChange', 1000000);

    [fr , gf]= fit(Time,dynamics',yfit1,opts)
    
    global coef_B2
    global coef_tau2
    
    coef_A1 = fr.A1;
    coef_A2 = fr.A2;
    coef_A3 = fr.A3;
    coef_B1 = fr.B1;
    coef_B2 = fr.B2;
    coef_c = fr.c;
    coef_s = fr.s;
    coef_tau1 = fr.tau1;
    coef_tau2 = fr.tau2;


    axes(handles.axes6);
    delete(h4);
    delete(h3);
    h4 = plot(handles.axes6, Time, coef_A1*exp(-(Time-coef_c).^2 /(2*coef_s^2)) + coef_A2*(coef_c-Time)/coef_s^3.*exp(-(Time-coef_c).^2/(2*coef_s^2)) + coef_A3*(Time.^2-2*coef_c*Time-coef_s^2 + coef_c^2)/coef_s^5.*exp(-(Time-coef_c).^2 /(2*coef_s^2)) +coef_B1*exp(-Time/coef_tau1).*exp( (coef_c+coef_s^2/(2*coef_tau1))/coef_tau1).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau1))/(1.41421356237*coef_s)))+coef_B2*exp(-Time/coef_tau2).*exp( (coef_c+coef_s^2/(2*coef_tau2))/coef_tau2).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau2))/(1.41421356237*coef_s))),'DisplayName','Fit result','color',[1 0 0]  );
    axis([t1 t2 -inf inf]);
    set([h4],'LineWidth',2);  
    grid on;
    
    msgbox( sprintf('ARTIFACT:\nA1 = %d\nA2 = %d\nA3 = %d\n\nEXPONENTIAL AMPLITUDES:\nB1 = %d\nB2 = %d\n\nc = %d\ns = %d\n\nDECAY CONSTANTS:\ntau_1 = %d\ntau_2 = %d',coef_A1,coef_A2,coef_A3,coef_B1,coef_B2,coef_c,coef_s,coef_tau1,coef_tau2));
    msgbox( sprintf('GOODNES OF FIT:\nrsquare = %.3d\nadjrsquare = %.3d\nrmse = %.3d',gf.rsquare,gf.adjrsquare,gf.rmse));
elseif pop ==3
    yfit1 = fittype('A1*exp(-(t-c).^2 /(2*s^2)) + A2*(c-t)/s^3.*exp(-(t-c).^2/(2*s^2)) + A3*(t.^2-2*c*t-s^2 + c^2)/s^5.*exp(-(t-c).^2 /(2*s^2)) +B1*exp(-t/tau1).*exp( (c+s^2/(2*tau1))/tau1).*(1+erf( (t-c-(s^2/tau1))/(1.41421356237*s)))+B2*exp(-t/tau2).*exp( (c+s^2/(2*tau2))/tau2).*(1+erf( (t-c-(s^2/tau2))/(1.41421356237*s)))+B3*exp(-t/tau3).*exp( (c+s^2/(2*tau3))/tau3).*(1+erf( (t-c-(s^2/tau3))/(1.41421356237*s)))' ,'ind','t'); 
    opts=fitoptions(yfit1);
    set(opts,'robust','off','Algorithm', 'Levenberg-Marquardt','tolfun',1e-18,'tolx',1e-18,'maxiter',10000,'MaxFunEvals',10000, ...
        'StartPoint',[A1 A2 A3 B1 B2 B3 c s tau1 tau2 tau3 ],'Lower', [A1_min A2_min A3_min B1_min B2_min B3_min c_min s_min tau_1_min tau_2_min tau_3_min],'Upper', [A1_max A2_max A3_max B1_max B2_max B3_max c_max s_max tau_1_max tau_2_max tau_3_max], 'DiffMinChange', 1.0000e-08, 'DiffMaxChange', 1000000);

    [fr , gf]= fit(Time,dynamics',yfit1,opts)
    
    global coef_B2
    global coef_tau2
    global coef_B3
    global coef_tau3
    
    coef_A1 = fr.A1;
    coef_A2 = fr.A2;
    coef_A3 = fr.A3;
    coef_B1 = fr.B1;
    coef_B2 = fr.B2;
    coef_B3 = fr.B3;
    coef_c = fr.c;
    coef_s = fr.s;
    coef_tau1 = fr.tau1;
    coef_tau2 = fr.tau2;
    coef_tau3 = fr.tau3;
    
    
    axes(handles.axes6);
    delete(h4);
    delete(h3);
    h4 = plot(handles.axes6, Time, coef_A1*exp(-(Time-coef_c).^2 /(2*coef_s^2)) + coef_A2*(coef_c-Time)/coef_s^3.*exp(-(Time-coef_c).^2/(2*coef_s^2)) + coef_A3*(Time.^2-2*coef_c*Time-coef_s^2 + coef_c^2)/coef_s^5.*exp(-(Time-coef_c).^2 /(2*coef_s^2)) +coef_B1*exp(-Time/coef_tau1).*exp( (coef_c+coef_s^2/(2*coef_tau1))/coef_tau1).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau1))/(1.41421356237*coef_s)))+coef_B2*exp(-Time/coef_tau2).*exp( (coef_c+coef_s^2/(2*coef_tau2))/coef_tau2).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau2))/(1.41421356237*coef_s)))+coef_B3*exp(-Time/coef_tau3).*exp( (coef_c+coef_s^2/(2*coef_tau3))/coef_tau3).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau3))/(1.41421356237*coef_s))),'DisplayName','Fit result','color',[1 0 0]  );
    axis([t1 t2 -inf inf]);
    set([h4],'LineWidth',2);  
    grid on;
    
    msgbox( sprintf('ARTIFACT:\nA1 = %d\nA2 = %d\nA3 = %d\n\nEXPONENTIAL AMPLITUDES:\nB1 = %d\nB2 = %d\nB3 = %d\n\nc = %d\ns = %d\n\nDECAY CONSTANTS:\ntau_1 = %d\ntau_2 = %d\ntau_3 = %d',coef_A1,coef_A2,coef_A3,coef_B1,coef_B2,coef_B3,coef_c,coef_s,coef_tau1,coef_tau2,coef_tau3));
    msgbox( sprintf('GOODNES OF FIT:\nrsquare = %.3d\nadjrsquare = %.3d\nrmse = %.3d',gf.rsquare,gf.adjrsquare,gf.rmse));
elseif pop ==4
    yfit1 = fittype('A1*exp(-(t-c).^2 /(2*s^2)) + A2*(c-t)/s^3.*exp(-(t-c).^2/(2*s^2)) + A3*(t.^2-2*c*t-s^2 + c^2)/s^5.*exp(-(t-c).^2 /(2*s^2)) +B1*exp(-t/tau1).*exp( (c+s^2/(2*tau1))/tau1).*(1+erf( (t-c-(s^2/tau1))/(1.41421356237*s)))+B2*exp(-t/tau2).*exp( (c+s^2/(2*tau2))/tau2).*(1+erf( (t-c-(s^2/tau2))/(1.41421356237*s)))+B3*exp(-t/tau3).*exp( (c+s^2/(2*tau3))/tau3).*(1+erf( (t-c-(s^2/tau3))/(1.41421356237*s)))+B4*exp(-t/tau4).*exp( (c+s^2/(2*tau4))/tau4).*(1+erf( (t-c-(s^2/tau4))/(1.41421356237*s)))' ,'ind','t'); 
    opts=fitoptions(yfit1);
    set(opts,'robust','off','Algorithm', 'Levenberg-Marquardt','tolfun',1e-18,'tolx',1e-18,'maxiter',10000,'MaxFunEvals',10000, ...
        'StartPoint',[A1 A2 A3 B1 B2 B3 B4 c s tau1 tau2 tau3 tau4],'Lower', [A1_min A2_min A3_min B1_min B2_min B3_min B4_min c_min s_min tau_1_min tau_2_min tau_3_min tau_4_min],'Upper', [A1_max A2_max A3_max B1_max B2_max B3_max B4_max c_max s_max tau_1_max tau_2_max tau_3_max tau_4_max], 'DiffMinChange', 1.0000e-08, 'DiffMaxChange', 1000000);

    [fr , gf]= fit(Time,dynamics',yfit1,opts)
    
    
    global coef_B2
    global coef_tau2
    global coef_B3
    global coef_tau3
    global coef_B4
    global coef_tau4
    

    coef_A1 = fr.A1;
    coef_A2 = fr.A2;
    coef_A3 = fr.A3;
    coef_B1 = fr.B1;
    coef_B2 = fr.B2;
    coef_B3 = fr.B3;
    coef_B4 = fr.B4;
    coef_c = fr.c;
    coef_s = fr.s;
    coef_tau1 = fr.tau1;
    coef_tau2 = fr.tau2;
    coef_tau3 = fr.tau3;
    coef_tau4 = fr.tau4;
    
    axes(handles.axes6);
    delete(h4);
    delete(h3);
    h4 = plot(handles.axes6, Time, coef_A1*exp(-(Time-coef_c).^2 /(2*coef_s^2)) + coef_A2*(coef_c-Time)/coef_s^3.*exp(-(Time-coef_c).^2/(2*coef_s^2)) + coef_A3*(Time.^2-2*coef_c*Time-coef_s^2 + coef_c^2)/coef_s^5.*exp(-(Time-coef_c).^2 /(2*coef_s^2)) +coef_B1*exp(-Time/coef_tau1).*exp( (coef_c+coef_s^2/(2*coef_tau1))/coef_tau1).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau1))/(1.41421356237*coef_s)))+coef_B2*exp(-Time/coef_tau2).*exp( (coef_c+coef_s^2/(2*coef_tau2))/coef_tau2).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau2))/(1.41421356237*coef_s)))+coef_B3*exp(-Time/coef_tau3).*exp( (coef_c+coef_s^2/(2*coef_tau3))/coef_tau3).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau3))/(1.41421356237*coef_s)))+coef_B4*exp(-Time/coef_tau4).*exp( (coef_c+coef_s^2/(2*coef_tau4))/coef_tau4).*(1+erf( (Time-coef_c-(coef_s^2/coef_tau4))/(1.41421356237*coef_s))),'DisplayName','Fit result','color',[1 0 0]  );
    axis([t1 t2 -inf inf]);
    set([h4],'LineWidth',2);  
    grid on;

    msgbox( sprintf('ARTIFACT:\nA1 = %d\nA2 = %d\nA3 = %d\n\nEXPONENTIAL AMPLITUDES:\nB1 = %d\nB2 = %d\nB3 = %d\nB4 = %d\n\nc = %d\ns = %d\n\nDECAY CONSTANTS:\ntau_1 = %d\ntau_2 = %d\ntau_3 = %d\ntau_4 = %d',coef_A1,coef_A2,coef_A3,coef_B1,coef_B2,coef_B3,coef_B4,coef_c,coef_s,coef_tau1,coef_tau2,coef_tau3,coef_tau4));
    msgbox( sprintf('GOODNES OF FIT:\nrsquare = %.3d\nadjrsquare = %.3d\nrmse = %.3d',gf.rsquare,gf.adjrsquare,gf.rmse));
end

fit_line = fr(Time);
XPM = dynamics*0;
XPM = XPM';
a = size(Time);
for i = 1:a(1)
    XPM(i) = coef_A1*exp(-(Time(i)-coef_c)^2 /(2*coef_s^2)) + coef_A2*(coef_c-Time(i))/coef_s^3* exp(-(Time(i)-coef_c)^2 /(2*coef_s^2))+ coef_A3*((Time(i))^2-2*coef_c*Time(i)-coef_s^2 + coef_c^2)/coef_s^5*exp(-(Time(i)-coef_c)^2 /(2*coef_s^2));
    i;
end


%--------------------------------------------------------------------------
%                           CALLBACKS
%--------------------------------------------------------------------------

function Time_1_Callback(~, ~, handles)
global a
global pop_up
global valna1
global checker_fit

if checker_fit == 0
    cla(handles.axes4);
    for r = 1:a
        plot_svd_dynamics(r,handles);
        hold on
    end
end
cla(handles.axes6);
plot_cut(valna1,handles);
hold on;
plot_ansatz(pop_up,handles)
plot_fit(pop_up,handles)

% plot2d(handles) %uncomment for 2d plot variable scale

% --- Executes during object creation, after setting all properties.
function Time_1_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Time_2_Callback(~, ~, handles)
global a
global pop_up
global valna1
global checker_fit

if checker_fit == 0
    cla(handles.axes4);
    for r = 1:a
        plot_svd_dynamics(r,handles);
        hold on
    end
end
cla(handles.axes6);
plot_cut(valna1,handles);
hold on;
plot_ansatz(pop_up,handles)
plot_fit(pop_up,handles)

% plot2d(handles) %uncomment for 2d plot variable scale

% --- Executes during object creation, after setting all properties.
function Time_2_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Wavelength1_Callback(~, ~, handles)
global a
cla(handles.axes5);
for r = 1:a
    plot_svd_spectrum(r,handles);
    hold on
end

% cla(handles.axes6);
% plot_cut(valna1,handles);
% hold on;
% plot_ansatz(pop_up,handles)
% plot_fit(pop_up,handles)

% plot2d(handles) %odkomentirati za mijenjanje skale 2d matrice

% --- Executes during object creation, after setting all properties.
function Wavelength1_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Wavelength2_Callback(~, ~, handles)
global a
cla(handles.axes5);
for r = 1:a
    plot_svd_spectrum(r,handles);
    hold on
end

% cla(handles.axes6);
% plot_cut(valna1,handles);
% hold on;
% plot_ansatz(pop_up,handles)
% plot_fit(pop_up,handles)

% plot2d(handles) %odkomentirati za mijenjanje skale 2d matrice

function Wavelength2_CreateFcn(hObject, ~, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function C1_Callback(~, ~, handles)
plot2d(handles)

function C1_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function C2_Callback(~, ~, handles)
plot2d(handles)
function C2_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SVD_value_number_Callback(~, ~, handles)

global a
a = str2double(get(handles.SVD_value_number,'String'));

%plotting dynamics
cla(handles.axes4);
cla(handles.axes5);
for r = 1:a
    plot_svd_dynamics(r,handles);
    hold on
    plot_svd_spectrum(r,handles);
    hold on
end



function SVD_value_number_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function axes3_CreateFcn(~, ~, ~)

function Wavelength_fit_Callback(~, ~, handles)
global valna1

valna1 = str2double(get(handles.Wavelength_fit,'String'));
plot_cut(valna1,handles);
hold on;

function Wavelength_fit_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(~, ~, handles)
global pop_up
pop_up = get(handles.popupmenu3,'Value');
plot_ansatz(pop_up,handles)

function popupmenu3_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function B1_Callback(~, ~, handles)
global B1
B1 = str2double(get(handles.B1,'String'));
global pop_up

pop_up = get(handles.popupmenu3,'Value');
plot_ansatz(pop_up,handles)

% --- Executes during object creation, after setting all properties.
function B1_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function c_Callback(~, ~, handles)
global c
c = str2double(get(handles.c,'String'))
global pop_up
pop_up = get(handles.popupmenu3,'Value');
plot_ansatz(pop_up,handles)

function c_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s_Callback(~, ~, handles)
global s
s = str2double(get(handles.s,'String'))
global pop_up
pop_up = get(handles.popupmenu3,'Value');
plot_ansatz(pop_up,handles)

function s_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tau1_Callback(~, ~, handles)
global tau1
tau1 = str2double(get(handles.tau1,'String'));
global pop_up

pop_up = get(handles.popupmenu3,'Value');
plot_ansatz(pop_up,handles)

function tau1_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function B2_Callback(~, ~, handles)

global B2
B2 = str2double(get(handles.B2,'String'));
global pop_up

pop_up = get(handles.popupmenu3,'Value');
plot_ansatz(pop_up,handles)

function B2_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tau_2_Callback(~, ~, handles)
global tau2
tau2 = str2double(get(handles.tau_2,'String'));
global pop_up
pop_up = get(handles.popupmenu3,'Value');
plot_ansatz(pop_up,handles)

function tau_2_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function B3_Callback(~, ~, handles)

global B3
B3 = str2double(get(handles.B3,'String'));
global pop_up

pop_up = get(handles.popupmenu3,'Value');
plot_ansatz(pop_up,handles)

function B3_CreateFcn(hObject, ~, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function B4_Callback(~, ~, handles)
global B4
B4 = str2double(get(handles.B4,'String'));
global pop_up

pop_up = get(handles.popupmenu3,'Value');
plot_ansatz(pop_up,handles)

function B4_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Tau_3_Callback(~, ~, handles)
global tau3
tau3 = str2double(get(handles.Tau_3,'String'));
global pop_up

pop_up = get(handles.popupmenu3,'Value');
plot_ansatz(pop_up,handles)

function Tau_3_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Tau_4_Callback(~, ~, handles)
global tau4
tau4 = str2double(get(handles.Tau_4,'String'));
global pop_up

pop_up = get(handles.popupmenu3,'Value');
plot_ansatz(pop_up,handles)

function Tau_4_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function A1_Callback(~, ~, handles)

global A1
A1 = str2double(get(handles.A1,'String'));
global pop_up

pop_up = get(handles.popupmenu3,'Value');
plot_ansatz(pop_up,handles)

function A1_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function A2_Callback(~, ~, handles)
global A2
A2 = str2double(get(handles.A2,'String'));
global pop_up

pop_up = get(handles.popupmenu3,'Value');
plot_ansatz(pop_up,handles)

function A2_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function A3_Callback(~, ~, handles)
global A3

A3 = str2double(get(handles.A3,'String'));
global pop_up

pop_up = get(handles.popupmenu3,'Value');
plot_ansatz(pop_up,handles)

function A3_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FIT CURVE.
function fit_button_Callback(~, ~, handles)

pop_up = get(handles.popupmenu3,'Value');

fit_and_plot(pop_up, handles)


function c_min_Callback(~, ~, handles)
global c_min
c_min = str2double(get(handles.c_min,'String'));

function c_min_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function c_max_Callback(~, ~, handles)
global c_max
c_max = str2double(get(handles.c_max,'String'));

function c_max_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s_min_Callback(~, ~, handles)
global s_min
s_min = str2double(get(handles.s_min,'String'));

function s_min_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s_max_Callback(~, ~, handles)
global s_max
s_max = str2double(get(handles.s_max,'String'));

function s_max_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function A1_min_Callback(~, ~, handles)
global A1_min
A1_min = str2double(get(handles.A1_min,'String'))

function A1_min_CreateFcn(hObject, ~, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function A2_min_Callback(~, ~, handles)
global A2_min
A2_min = str2double(get(handles.A2_min,'String'));

function A2_min_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function A3_min_Callback(~, ~, handles)
global A3_min
A3_min = str2double(get(handles.A3_min,'String'));

function A3_min_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function A1_max_Callback(~, ~, handles)
global A1_max
A1_max = str2double(get(handles.A1_max,'String'));

function A1_max_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function A2_max_Callback(~, ~, handles)
global A2_max
A2_max = str2double(get(handles.A2_max,'String'));

function A2_max_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function A3_max_Callback(~, ~, handles)
global A3_max
A3_max = str2double(get(handles.A3_max,'String'));
function A3_max_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function B1_min_Callback(~, ~, handles)
global B1_min
B1_min = str2double(get(handles.B1_min,'String'));

function B1_min_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function B2_min_Callback(~, ~, handles)
global B2_min
B2_min = str2double(get(handles.B2_min,'String'));

function B2_min_CreateFcn(hObject, ~, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function B3_min_Callback(~, ~, handles)
global B3_min
B3_min = str2double(get(handles.B3_min,'String'));

function B3_min_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function B4_min_Callback(~, ~, handles)
global B4_min
B4_min = str2double(get(handles.B4_min,'String'));

function B4_min_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function B1_max_Callback(~, ~, handles)
global B1_max
B1_max = str2double(get(handles.B1_max,'String'));

function B1_max_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function B2_max_Callback(~, ~, handles)
global B2_max
B2_max = str2double(get(handles.B2_max,'String'));

function B2_max_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function B3_max_Callback(~, ~, handles)
global B3_max
B3_max = str2double(get(handles.B3_max,'String'));

function B3_max_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function B4_max_Callback(~, ~, handles)
global B4_max
B4_max = str2double(get(handles.B4_max,'String'));

function B4_max_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tau_1_min_Callback(~, ~, handles)
global tau_1_min
tau_1_min = str2double(get(handles.tau_1_min,'String'));

function tau_1_min_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tau_2_min_Callback(~, ~, handles)
global tau_2_min
tau_2_min = str2double(get(handles.tau_2_min,'String'));

function tau_2_min_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tau_3_min_Callback(~, ~, handles)
global tau_3_min
tau_3_min = str2double(get(handles.tau_3_min,'String'));

function tau_3_min_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tau_4_min_Callback(~, ~, handles)
global tau_4_min
tau_4_min = str2double(get(handles.tau_4_min,'String'));

function tau_4_min_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tau_1_max_Callback(~, ~, handles)
global tau_1_max
tau_1_max = str2double(get(handles.tau_1_max,'String'));

function tau_1_max_CreateFcn(hObject, ~, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tau_2_max_Callback(~, ~, handles)
global tau_2_max
tau_2_max = str2double(get(handles.tau_2_max,'String'));

function tau_2_max_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tau_3_max_Callback(~, ~, handles)
global tau_3_max
tau_3_max = str2double(get(handles.tau_3_max,'String'));

function tau_3_max_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tau_4_max_Callback(~, ~, handles)
global tau_4_max
tau_4_max = str2double(get(handles.tau_4_max,'String'));

function tau_4_max_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in plot_residual.
function plot_residual_Callback(~, ~, handles)
global fit_line
global dynamics
global axesfont

figure(5)
size(dynamics)
size(fit_line)
Residual = dynamics' - fit_line;
h5 = plot(Residual);
set([h5],'LineWidth',2)
hold on;
grid on;
xlabel(['Delay [fs]'],'FontSize',axesfont);
ylabel('\Delta A','FontSize',axesfont);
set(gca, 'FontSize', 15);box on;grid on;
axis([-inf inf -inf inf]);


% --- Executes on button press in SUBTRACT ARTIFACT.
function subtracted_artifact_Callback(~, ~, handles)
global fit_line
global dynamics
global axesfont
global Time
global XPM


figure(6)
h6 = plot(Time,fit_line - XPM,'DisplayName','Fit - Artifact');
set([h6],'LineWidth',3)
hold on;
h7 = plot(Time,dynamics'-XPM,'DisplayName','Measurement - Artifact');
set([h7],'LineWidth',3)
xlabel(['Delay [fs]'],'FontSize',axesfont);
ylabel('\Delta A','FontSize',axesfont);
xt = get(gca, 'XTick');set(gca, 'FontSize', 15);box on;grid on;
grid on;
lgd = legend('Location','best');
axis([-inf 5000 -inf inf]);
