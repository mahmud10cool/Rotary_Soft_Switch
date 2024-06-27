clear; clc;
%% Outlet
% Parameters        
piston_dia = 8e-3;                        %
r = 20e-3;
Ap = 0.25 * pi * piston_dia^2;
alpha = deg2rad(45);                       % Swash plate angle
omega_rpm = 3600;
omega = (2*pi/60) * omega_rpm;
t = 0:0.0001:0.1;
theta_n = omega * t;
beta = 1.6e9;
P1 = 21e6;

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
        Q1(i) = 0;
    end
end

% Plot for the flow rate
figure(1)
plot(t, Q1)
title('Outlet Flow Rate vs. time')
xlabel('s')
ylabel('m^3/s')
