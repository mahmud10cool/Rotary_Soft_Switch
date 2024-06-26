% Author: Mahmud Suhaimi Ibrahim
% Rotary system dynamic programming code

%% Dynamic Programming for Rotary Soft Switch
clc; clear;
format long;

% State, Input and Time array sizes
P1_size = 11;
P2_size = 11;
thetadot_size = 201;
Av_size = 1;
Telec_size = 1;
FINALTIME = 0.15;
STEPSIZE = 1e-5;
TIMESIZE = 1 + FINALTIME/STEPSIZE;
weight = 5;
mm = 1; % number of minor steps

save ARRAY.mat P1_size P2_size thetadot_size Av_size Telec_size FINALTIME
arr_size = [num2str(P1_size), ' ', num2str(P2_size), ' ', ...
    num2str(thetadot_size), ' ', num2str(Av_size), ' ', num2str(Telec_size), ...
    ' ', num2str(FINALTIME)];

disp(arr_size)

% Maximum and minimum areas
Ap = 0.25*pi*(20e-3)^2;
Av_max = 0.5*Ap;

% Maximum and minimum torques
Tmin = 0;
Tmax = 0;

% Maximum and minimum angular velocities
theta_dot_min = -200;
theta_dot_max = 200;

%% Parameters
% Volume of chambers 
V1_0 = 1e-3;
V2_0 = 20e-6;

% Pressure rails
P_H = 20e6;
P_M = 10e6;

% Valve things
max_Avt = 0.5*0.25*pi*(20e-3)^2;
Cd = 0.6;

% Fluid properties
beta = 1.8e9;
rho = 870;

% Electric motor stuff
Kt = 70.5e-3;
Ke = 70.5e-3;
J_elec = 140e-7;

% Pump things
J_hyd = 3000e-7;
hyd_D = 1.6e-6;    % In cc/rev

%% Time discretization

tf = FINALTIME;
dt = tf/(TIMESIZE-1);
t = 0:dt:tf;
n_t = length(t);

%% State discretization

% Set the pressure resolution for the pressure arrays
n_P1 = P1_size;
P1 = linspace(0.75*P_M, 1.625*P_H, n_P1);

if ~any(P1==P_H) 
    return
elseif ~any(P1==P_M)
    return
end

n_P2 = P2_size;
P2 = linspace(0.75*P_M, 1.625*P_H, n_P2);

if ~any(P2==P_H) 
    return
elseif ~any(P2==P_M)
    return
end

% Speed discretization
n_theta_dot = thetadot_size;
theta_dot = linspace(theta_dot_min, theta_dot_max, n_theta_dot);

if ~any(theta_dot==0)
    return
end

[P1_matrix, P2_matrix, theta_dot_matrix] = ndgrid(P1,P2,theta_dot);

P1_states = P1_matrix(:);
P2_states = P2_matrix(:);
theta_dot_states = theta_dot_matrix(:);

% Number of discretizations by changing ndgrid matrix into column vector
States = length(P1_states);

%% Input discretization

n_Av = Av_size;
A_vt = linspace(0, Av_max, n_Av);

n_Telec = Telec_size;
Telec = linspace(Tmin, Tmax, n_Telec);

[A_vt_matrix, Telec_matrix] = ndgrid(A_vt, Telec);

A_vt_options = A_vt_matrix(:)';
Telec_options = Telec_matrix(:)';

U_input = length(A_vt_options);

%% Dynamic Programming

% Initialising the cost matrix
J_cost = NaN(n_t,States);

% This is the thing altering the results
ind_final_P1 = find(P1 == P_H);
ind_final_P2 = find(P2 == P_H);
ind_final_theta_dot = find(theta_dot == 0);

ind_final = sub2ind(size(P1_matrix), ind_final_P1, ind_final_P2, ind_final_theta_dot);

%  Cost at final time
J_cost(end,:) = weight*(1e-5*abs(P1_states-P1(ind_final_P1)).^1 + 1e-5*abs(P2_states - P2(ind_final_P2)).^1 + ...
    10*abs((theta_dot_states-theta_dot(ind_final_theta_dot))).^1);
% J_cost(end,:) = 0;

% Initialising the optimal input (index) matrix
u_opt_matrix = NaN(n_t,States); 
u_opt_matrix(end,:) = 1;


%% Precalculations

% I will use 'k' for P2, 'l' for P1 and 'm' for theta_dot

%% First term
Q1 = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-P2_states)).*sign(P_H-P2_states);

% Current revolutions
theta_1 = ((P1_states-P_M)*V1_0*2*pi)/(beta*hyd_D);

% Flow to the chamber via the pump/motor
Q_hyd_1 = (hyd_D/2*pi)*theta_dot_states;

% Derivative of chamber pressure
P1dot = (beta./ V1_0) .* (Q_hyd_1);

% Derivative of small pressure
P2dot = (beta ./ V2_0) .* (Q1 - Q_hyd_1);

% Acceleration calculations
hyd_torque_1 = (hyd_D/(2*pi))*(P2_states - P1_states);
theta_doubledot = (1/(J_hyd + J_elec))*(hyd_torque_1 - Telec_options);

k1 = P2dot;
l1 = P1dot;
m1 = theta_doubledot;

%% Second term
Q2 = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-(P2_states+0.5*dt*k1))).*sign(P_H-(P2_states+0.5*dt*k1));

% Current revolutions
theta_2 = (((P1_states+0.5*dt*l1)-P_M)*V1_0*2*pi)/(beta*hyd_D);

% Flow to the chamber via the pump/motor
Q_hyd_2 = (hyd_D/2*pi)*(theta_dot_states+0.5*dt*m1);

% Hydraulic Torque
hyd_torque_2 = (hyd_D/(2*pi))*((P2_states+0.5*dt*k1) - (P1_states+0.5*dt*l1));

k2 = (beta ./ V2_0) .* (Q2 - Q_hyd_2);
l2 = (beta./ V1_0) .* (Q_hyd_2);
m2 = (1/(J_hyd + J_elec))*(hyd_torque_2 - Telec_options);

%% Third term
Q3 = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-(P2_states+0.5*dt*k2))).*sign(P_H-(P2_states+0.5*dt*k2));

% Current revolutions
theta_3 = (((P1_states+0.5*dt*l2)-P_M)*V1_0*2*pi)/(beta*hyd_D);

% Flow to the chamber via the pump/motor
Q_hyd_3 = (hyd_D/2*pi)*(theta_dot_states+0.5*dt*m2);

% Hydraulic Torque
hyd_torque_3 = (hyd_D/(2*pi))*((P2_states+0.5*dt*k2) - (P1_states+0.5*dt*l2));

k3 = (beta ./ V2_0) .* (Q3 - Q_hyd_3);
l3 = (beta./ V1_0) .* (Q_hyd_3);
m3 = (1/(J_hyd + J_elec))*(hyd_torque_3 - Telec_options);

%% Fourth term
Q4 = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-(P2_states+dt*k3))).*sign(P_H-(P2_states+dt*k3));

% Current revolutions
theta_4 = (((P1_states+dt*l3)-P_M)*V1_0*2*pi)/(beta*hyd_D);

% Flow to the chamber via the pump/motor
Q_hyd_4 = (hyd_D/2*pi)*(theta_dot_states+dt*m3);

% Hydraulic Torque
hyd_torque_4 = (hyd_D/(2*pi))*((P2_states+dt*k3) - (P1_states+dt*l3));

k4 = (beta ./ V2_0) .* (Q4 - Q_hyd_4);
l4 = (beta./ V1_0) .* (Q_hyd_4);
m4 = (1/(J_hyd + J_elec))*(hyd_torque_4 - Telec_options);

%% Integration
P2next = P2_states + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
P1next = P1_states + (dt/6)*(l1 + 2*l2 + 2*l3 + l4);
theta_dot_next = theta_dot_states + (dt/6)*(m1 + 2*m2 + 2*m3 + m4);

%% Throttling loss calculations
throttling_1 = Q1 .* (P_H-P2_states);

throttling_2 = Q2 .* (P_H-(P2_states+0.5*dt*k1));

throttling_3 = Q3 .* (P_H-(P2_states+0.5*dt*k2));

throttling_4 = Q4 .* (P_H-(P2_states+dt*k3));

energy_throttle_loss = (dt/6)*(throttling_1+2*throttling_2+2*throttling_3+throttling_4);

%% Starting the DP
startTimer = tic;

disp('DP started successfully.')

for t_ind = n_t-1:-1:1
    J = reshape(J_cost(t_ind+1,:), n_P1, n_P2, n_theta_dot);

    cost_next = interpn(P1, P2, theta_dot, J, P1next, P2next, theta_dot_next, 'linear');

    cost_next(isnan(cost_next)) = max(J(:));

    total_loss = cost_next + energy_throttle_loss;

    [J_cost(t_ind,:), u_opt_matrix(t_ind,:)] = min(total_loss, [], 2, 'omitnan');

   if rem((n_t-t_ind)/(n_t-1)*100,1) == 0
    disp([num2str(round((n_t-t_ind)/(n_t-1)*100)) '% Complete']);
   end

end

disp('DP completed successfully.')

%% Timing the DP

% Stopping the time
stopTimer = toc(startTimer);

% Finding out hours, minutes and seconds
timeTaken = seconds(stopTimer);
timeTaken.Format = 'hh:mm:ss';

% Displaying the time
printTime = ['Time taken to complete: ', string(timeTaken)];
disp(printTime)

%% Forward moving

% States vs time arrays
P2_f = NaN(n_t,1);
P1_f = NaN(n_t,1);
theta_dot_f = NaN(n_t,1);

% Flow arrays
Q_forward = NaN(n_t,1);
Q_next_forward = NaN(n_t,1);
Q_hyd_1_f = NaN(n_t,1);

% Find the right initial states
ind_P1 = find(P1 == P_M);
ind_P2 = find(P2 == P_M);
ind_theta_dot = find(theta_dot == 0);

% Initial states
P1_f(1,1) = P1(ind_P1);
P2_f(1,1) = P2(ind_P2);
theta_dot_f(1,1) = theta_dot(ind_theta_dot);

theta_1_f = NaN(n_t,1);
theta_2_f = NaN(n_t,1);
theta_3_f = NaN(n_t,1);
theta_4_f = NaN(n_t,1);

J_cost_f = NaN(n_t, 1);
u_opt_f = NaN(n_t, 1);
loss_f = NaN(n_t, 1);
total_loss_f = NaN(n_t,1);
total_loss_f(1) = 0;

% Input vs time arrays
Telec_f = NaN(n_t,1);
A_vt_f = NaN(n_t,1);

for t_f = 1:1:n_t-1

    J_f = reshape(J_cost(t_f+1,:), n_P1, n_P2, n_theta_dot);

    %% First term
    Q1_f = Cd*A_vt_options*sqrt((2/rho)*abs(P_H-P2_f(t_f)))*sign(P_H-P2_f(t_f));

    Q_hyd_1_f(t_f) = (hyd_D/(2*pi))*theta_dot_f(t_f);

    hyd_torque_1_f = (hyd_D/(2*pi))*(P2_f(t_f) - P1_f(t_f));
    
    P2dot_f = (beta/(V2_0))*(Q1_f - Q_hyd_1_f(t_f));
    P1dot_f = (beta/(V1_0)) * (Q_hyd_1_f(t_f));
    theta_doubledot_f = (1/(J_hyd + J_elec))*(hyd_torque_1_f - Telec_options);

    k1_f = P2dot_f;
    l1_f = P1dot_f;
    m1_f = theta_doubledot_f;

    %% Second term
    Q2_f = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-(P2_f(t_f)+0.5*dt*k1_f))).*sign(P_H-(P2_f(t_f)+0.5*dt*k1_f));

    Q_hyd_2_f = (hyd_D/(2*pi))*(theta_dot_f(t_f)+0.5*dt*m1_f);

    hyd_torque_2_f = (hyd_D/(2*pi))*((P2_f(t_f)+0.5*dt*k1_f) - (P1_f(t_f)+0.5*dt*l1_f));

    k2_f = (beta/(V2_0))*(Q2_f - Q_hyd_2_f);
    l2_f = (beta/(V1_0)) * (Q_hyd_2_f);
    m2_f = (1/(J_hyd + J_elec))*(hyd_torque_2_f - Telec_options);

    %% Third term
    Q3_f = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-(P2_f(t_f)+0.5*dt*k2_f))).*sign(P_H-(P2_f(t_f)+0.5*dt*k2_f));

    Q_hyd_3_f = (hyd_D/(2*pi))*(theta_dot_f(t_f)+0.5*dt*m2_f);

    hyd_torque_3_f = (hyd_D/(2*pi))*((P2_f(t_f)+0.5*dt*k2_f) - (P1_f(t_f)+0.5*dt*l2_f));

    k3_f = (beta/(V2_0))*(Q3_f - Q_hyd_3_f);
    l3_f = (beta/(V1_0)) * (Q_hyd_3_f);
    m3_f = (1/(J_hyd + J_elec))*(hyd_torque_3_f - Telec_options);

    %% Fourth term
    Q4_f = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-(P2_f(t_f)+dt*k3_f))).*sign(P_H-(P2_f(t_f)+dt*k3_f));

    Q_hyd_4_f = (hyd_D/(2*pi))*(theta_dot_f(t_f)+dt*m3_f);

    hyd_torque_4_f = (hyd_D/(2*pi))*((P2_f(t_f)+dt*k3_f) - (P1_f(t_f)+dt*l3_f));

    k4_f = (beta/(V2_0))*(Q4_f - Q_hyd_4_f);
    l4_f = (beta/(V1_0)) * (Q_hyd_4_f);
    m4_f = (1/(J_hyd + J_elec))*(hyd_torque_4_f - Telec_options);

    %% Integration
    P2next_f = P2_f(t_f) + (dt/6)*(k1_f + 2*k2_f + 2*k3_f + k4_f);
    P1next_f = P1_f(t_f) + (dt/6)*(l1_f + 2*l2_f + 2*l3_f + l4_f);
    theta_dot_next_f = theta_dot_f(t_f) + (dt/6)*(m1_f + 2*m2_f + 2*m3_f + m4_f);

    %% Theta
    theta_next_f = ((P1next_f-P_M)*V1_0*2*pi)/(beta*hyd_D);

    %% Cost analysis
    cost_next_f = interpn(P1, P2, theta_dot, J_f, P1next_f,P2next_f,theta_dot_next_f,'linear');

    cost_next_f(isnan(cost_next_f)) = max(J(:));

    %% Throttling losses
    throttling_1_f = Q1_f .* (P_H-P2_f(t_f));

    throttling_2_f = Q2_f .* (P_H-(P2_f(t_f)+0.5*dt*k1_f));

    throttling_3_f = Q3_f .* (P_H-(P2_f(t_f)+0.5*dt*k2_f));

    throttling_4_f = Q4_f .* (P_H-(P2_f(t_f)+dt*k3_f));

    energy_throttle_loss_f = (dt/6)*(throttling_1_f+2*throttling_2_f+2*throttling_3_f+throttling_4_f);

    total_cost_f = cost_next_f + energy_throttle_loss_f;

    % Cost minimization
    % The heart of the forward path planning
    [J_cost_f(t_f),u_opt_f(t_f)] = min(total_cost_f, [], 'omitnan');

    % Setting up the optimal path
    P2_f(t_f+1) = P2next_f(u_opt_f(t_f));
    P1_f(t_f+1) = P1next_f(u_opt_f(t_f));

    Q_forward(t_f) = Q1_f(u_opt_f(t_f));
    Q_next_forward(t_f) = Q4_f(u_opt_f(t_f));

    theta_dot_f(t_f+1) = theta_dot_next_f(u_opt_f(t_f));
    theta_1_f(t_f+1) = theta_next_f(u_opt_f(t_f));

    A_vt_f(t_f) = A_vt_options(u_opt_f(t_f));
    Telec_f(t_f) = Telec_options(u_opt_f(t_f));

    loss_f(t_f) = energy_throttle_loss_f(u_opt_f(t_f));
    total_loss_f(t_f+1) = total_loss_f(t_f) + loss_f(t_f);

end

%% Plots
% General plots
plot_fig(1) = figure;

subplot(3,2,1)
plot(t, P2_f, 'LineWidth',2);
xlabel('Time (s)')
ylabel('Pressure (Pa)')
grid
title('Small Pressure')

subplot(3,2,2)
plot(t, P1_f, 'LineWidth',2);
xlabel('Time (s)')
ylabel('Pressure (Pa)')
grid
title('Chamber Pressure')

subplot(3,2,3)
plot(t, Telec_f, 'LineWidth',2);
xlabel('Time (s)')
ylabel('Force (N)')
grid
title('Electric Torque Applied')

subplot(3,2,4)
plot(t, A_vt_f, 'LineWidth',2);
xlabel('Time (s)')
ylabel('Area (m^2)')
grid
title('Valve Area Opening')
ylim([0 1.5*max(A_vt_f)])

subplot(3,2,5)
plot(t, theta_dot_f, 'LineWidth',2);
xlabel('Time (s)')
ylabel('ms^-1')
grid
title('Velocity')

subplot(3,2,6)
plot(t, theta_1_f, 'LineWidth',2);
xlabel('Time (s)')
ylabel('m')
grid
title('Angle')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

sgtitle('Runge Kutta Method')

%% Energy plot
plot_fig(2) = figure;

subplot(3,2,1)
plot(t, total_loss_f, 'LineWidth', 3);
xlabel('Time (s)')
ylabel('Energy (J)')
grid
title('Energy Loss')

subplot(3,2,2)
regen = cumsum(Telec_f.*theta_dot_f*dt);
plot(t, regen, 'LineWidth', 3);
xlabel('Time (s)')
ylabel('Energy (J)')
grid
title('Regenerated')

subplot(3,2,3)
energy_in = cumsum(P_H.*Q_forward*dt);
plot(t, energy_in, 'LineWidth', 3);
xlabel('Time (s)')
ylabel('Energy (J)')
grid
title('Energy In')

subplot(3,2,4)
work_out = cumsum(P1_f.*Q_hyd_1_f*dt);
plot(t, work_out, 'LineWidth', 3);
xlabel('Time (s)')
ylabel('Energy (J)')
grid
title('Work Out')

subplot(3,2,5)
compress_energy = cumsum(P2_f.*(Q_forward-Q_hyd_1_f)*dt);
plot(t, compress_energy, 'LineWidth', 3);
xlabel('Time (s)')
ylabel('Energy (J)')
grid
title('Compressed Energy')

subplot(3,2,6)
kinetic_energy = 0.5*(J_elec+J_hyd)*theta_dot_f.^2;
energy_out = total_loss_f + regen + work_out + compress_energy + kinetic_energy;
plot(t, energy_out, 'LineWidth', 3);
xlabel('Time (s)')
ylabel('Energy (J)')
grid
title('Energy Out')

sgtitle('Runge Kutta Method')

savefig(plot_fig,'Rotary_dp_rk_normal')
