clear; clc;
%% Outlet
% Parameters        
piston_dia = 8e-3;                      % Diameter of each piston
r = 20e-3;                              % Piston pitch radius 
Ap = 0.25 * pi * piston_dia^2;          % Area of each piston
alpha = deg2rad(45);                    % Swash plate angle
omega_rpm = 3600;                       % Rotational speed of pump
omega = (2*pi/60) * omega_rpm;          % Angular speed of pump
t = 0:10e-5:0.1;                        % Time
theta_n = omega * t;                    % Variation in angle
beta = 1.6e9;                           % Bulk Modulus
P1 = 21e6;                              % Outlet pressure

% Calculated parameters
n = size(theta_n);
Dp = zeros(n);
Q_kinematic = zeros(n);
Q_compress = zeros(n);
Q_leakage = zeros(n);
Q1 = zeros(n);

for i = 1:n(2)
    Dp(i) = Ap * r * tan(alpha) * sin(theta_n(i));
    % Kinematic
    Q_kinematic(i) = Dp(i) * omega;

    % Compressibility
    Q_compress(i) = Q_kinematic(i) * (P1/beta);

    % Total outlet/discharge flow rate
    Q1(i) = Q_kinematic(i) - Q_compress(i) - Q_leakage(i);
    if Q1(i) >=0
        Q1(i) = Q1(i);
    else
        Q1(i) = NaN;
    end
end

%% Inlet
% Parameters
P2 = 101e3;                               % Tank or inlet pressure

% Calculated parameters
Dp = zeros(n);
Q_kinematic = zeros(n);
Q_compress_2 = zeros(n);
Q2 = zeros(n);

for i = 1:n(2)
    Dp(i) = Ap * r * tan(alpha) * sin(theta_n(i));
    % Kinematic
    Q_kinematic(i) = Dp(i) * omega;

    % Compressibility
    Q_compress_2(i) = Q_kinematic(i) * (P2/beta);

    % Total outlet/discharge flow rate
    Q2(i) = Q_kinematic(i) + Q_compress_2(i);
    if Q2(i) <= 0
        Q2(i) = Q2(i);
    else
        Q2(i) = NaN;
    end
end

% Overall flow rate graph
figure(2)
plot(t,Q1,t,Q2)
title('Outlet Flow Rate vs. time')
xlabel('s')
ylabel('m^3/s')
legend('Qout', 'Qin')

