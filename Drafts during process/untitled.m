%% TUSAŞ INTERN PROJECT  @ 29/01/2025

% ALL UNITS ARE SI.

clc; clear; close all

%% General Parameters
g   = 9.81;
rho = 1.225;

%% Firstly, we need to initiate some design parameters.
W0 = 100;   % Gross (Maximum) Take-off Weight
Nb = 2;     % Number of blades

T = W0*g;   % Thrust

%% Main Rotor
D = 0.977 * W0^(0.308);               % Diameter
R = D/2;                              % Radius
A = pi*R^2;                           % Disc Area
chord = 0.0619*(W0^0.226)/(Nb^0.426);

Zcg = -0.5;                   % Z direction distance of rotor hub from CG
Xcg = 0;                      % X direction distance of rotor hub from CG
Ycg = 0;                      % Y direction distance of rotor hub from CG

R_twistd = 9;                 % Blade root angle degree
R_twist  = 9*(pi/180);        % Blade root angle radian

%% Tail Rotor
Dtr = 0.0886 * W0^(0.393);  % Diameter
Rtr = Dtr/2;                % Radius
Atr = pi*Rtr^2;             % Disc Area

%% Disc Loading
DL = 2.12 * (W0^(1/3) - 0.57) * g;  % From regression
%DL = W0*g/(pi*R^2)                 % Momentum theory

%% Total Power
Pto = 0.2928 * W0^(0.9043); % Total power at GTOW

%% Power Loading and Figure of Merit
% There are two equations for figure of merit, one is regression from 
% article, one is from momentum theory. When we arrange these together,
% figure of merit becomes:
FM = 10^(log10(T/(Pto*1000)) + 0.5*log10(DL)) / sqrt(2*rho);

%% Airframe
Lf = 0.824*D^(1.056);   % Fuselage length
Lrt = 1.09*D^(1.03);    % Overall length

%% Weight Prediction
We = 0.59*W0;
Wpl = 0.31*W0;

%% Maximum Speed
Vmax = 39.8*W0^(0.242);

%% Angular Velocity

Omega = (57.967*D^(0.6149))/R;

%% Rate of Climb
Vc = 142*W0^(0.157);

%% Engine (EMRAX 188)
Psea = 37; % kW

%% Results
data = {
    'INPUTS:', sprintf('');
    'MTOW (kg)', sprintf('%.0f', W0);
    'Blade Number', sprintf('%.0f', Nb);
    '', sprintf('');
    'OUTPUTS:', sprintf('');
    'Empty Weight (kg)', sprintf('%.1f', We);
    'Payload (kg)', sprintf('%.1f', Wpl);
    'Main Rotor Radius (m)', sprintf('%.2f', R);
    'Tail Rotor Radius (m)', sprintf('%.2f', Rtr);
    'Takeoff Total Power (kW)', sprintf('%.2f', Pto);
    'Fuselage Length (m)', sprintf('%.2f', Lf);
    'Airframe Length (m)', sprintf('%.2f', Lrt)
};

f = figure('Name', 'Helicopter Parameter Estimation',...
    'Position', [800, 400, 422, 400]);
uit = uitable(f, 'Data', data, ...
    'ColumnName', {'Parameter', 'Value'}, ...
    'RowName', [], ...
    'Units', 'normalized', ...
    'ColumnWidth', {320, 100}, ... % Sütun genişlikleri
    'Position', [0, 0, 1, 1],"FontSize",16);

%%
run("HelicopterParameterEstimation.mlapp");

% Vtip = (57.967*D^(0.6149));

% Vi = sqrt(T/2/rho/A);

% 
% Pi = T * Vi;
% Cpi = Pi/(rho*A*Vtip^3);
% 
% solid = Nb*c/pi/R;
% Cd0 = 0.008;


%%

% Mehmet Tunahan Kara ID:622942
% Celal Aksoy         ID:801550 


% 1- DESIGN TRENDS FOR ROTARY-WING UNMANNED AIR VEHICLES.
%    Khromov, V., & Rand, O.