clear; clc;
%% Outlet
% Parameters        
d_cyl = 8e-3;                                   % Diameter of each piston
r_p = 20e-3;                                    % Piston pitch radius 
A_cyl = 0.25 * pi * d_cyl^2;                    % Area of each piston
omega_rpm = 3600;                               % Rotational speed of pump
omega = (2*pi/60) * omega_rpm;                  % Angular speed of pump
beta = 1.6e9;                                   % Bulk Modulus
P1 = 21e6;                                      % Outlet pressure
S = 15e-3;                                      % Stroke length
Cd = 0.7;                                       % Coefficient of discharge
rho = 870;                                      % Density
P_T = 101325;                                   % Tank pressure
P_P = 21e6;                                     % System pressure

% Time parameters
theta_cycle = 2 * pi;
t_cycle = theta_cycle / omega;
theta_t = 2 * asin(d_cyl/(4*r_p));
t_t = theta_t/omega;
theta_d = theta_cycle/2 - 2 * theta_t; 
t_d = theta_d/omega;

% Valve area
step_size = 1e-7;

t1 = 0:step_size:t_t/2-step_size;
A_T1 = zeros(size(t1)); 

t2 = t_t/2:step_size:t_t-step_size;
A_T2 = ((2*A_cyl)/t_t) * (t2 - 0.5*t_t);

t3 = t_t:step_size:t_d+t_t-step_size;
A_T3 = A_cyl*ones(size(t3));

t4 = t_d+t_t:step_size:t_d+1.5*t_t-step_size;
A_T4 = -((2*A_cyl)/t_t) * (t4 - (t_d + 1.5*t_t));

t5 = t_d+1.5*t_t-step_size:step_size:t_cycle;
A_T5 = zeros(size(t5));

t = [t1 t2 t3 t4 t5];
A_T = [A_T1 A_T2 A_T3 A_T4 A_T5];

% Plot of valve area
figure(1)
plot(t, A_T, LineWidth=2, Color='blue')
title('Valve Port Area')
xlim([0, t_cycle])
xlabel('t / s')
ylabel('Area / m^2')

% Calculated parameters
Vdisp = 0.25 * pi * d_cyl^2 * S;                % Maximum volume displaced
Vo = Vdisp;                                     % Leftover volume
Vc = Vo + (Vdisp/2) * cos(omega*t) + Vdisp/2;   % The chamber volume
Vc_dot = -(Vdisp/2) * sin(omega*t);             % The change in the chamber volume.

% Plotting the volume change
figure(2)
plot(t, Vc, LineWidth=2, Color='blue')
title('Chamber Volume')
xlim([0, t_cycle])
xlabel('t / s')
ylabel('Vc / m^3')
