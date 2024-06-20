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
FINALTIME = 1e-3;
STEPSIZE = 1e-7;
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
V2_0 = 50e-6;

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
    disp('Stop 1');
    return
elseif ~any(P1==P_M)
    disp('Stop 2');
    return
end

n_P2 = P2_size;
P2 = linspace(0.75*P_M, 1.625*P_H, n_P2);

if ~any(P2==P_H)
    disp('Stop 3');
    return
elseif ~any(P2==P_M)
    disp('Stop 4');
    return
end

% Speed discretization
n_theta_dot = thetadot_size;
theta_dot = linspace(theta_dot_min, theta_dot_max, n_theta_dot);

if ~any(theta_dot==0)
    disp('Stop 5');
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
% Current revolutions
theta_now = ((P1_states-P_M)*V1_0*2*pi)/(beta*hyd_D);

% Flow to the chamber via the pump/motor
Q_hyd = (hyd_D/2*pi)*theta_dot_states;

% Chamber pressure
P1dot = (beta./ V1_0) .* (Q_hyd);
P1next = P1_states + P1dot * dt;
P1next = repmat(P1next, 1, U_input);

% Current flow rate
Qnow = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-P2_states)).*sign(P_H-P2_states);

% Derivative of small pressure
P2dot = (beta ./ V2_0) .* (Qnow - Q_hyd);
P2next = P2_states + P2dot*dt;

% Acceleration calculations
hyd_torque = (hyd_D/(2*pi))*(P2_states - P1_states);
theta_doubledot = (1/(J_hyd + J_elec))*(hyd_torque - Telec_options);
theta_dot_next = theta_dot_states + theta_doubledot * dt;

% Throttling loss calculations
throttling_now = Qnow.*(P_H-P2_states);
Qnext = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H - P2next)).*sign(P_H-P2next);
throttling_next = Qnext.*(P_H-P2next);
energy_throttle_loss = 0.5*(throttling_now + throttling_next)*dt;

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

% Find the right initial states
ind_P1 = find(P1 == P_M);
ind_P2 = find(P2 == P_M);
ind_theta_dot = find(theta_dot == 0);

% Initial states
P1_f(1,1) = P1(ind_P1);
P2_f(1,1) = P2(ind_P2);
theta_dot_f(1,1) = theta_dot(ind_theta_dot);

theta_f = NaN(n_t,1);
theta_f(1,1) = ((P1_f(1,1)-P_M)*V1_0*2*pi)/(beta*hyd_D);

J_cost_f = NaN(n_t, 1);
u_opt_f = NaN(n_t, 1);
loss_f = NaN(n_t, 1);
total_loss_f = NaN(n_t,1);
total_loss_f(1) = 0;

% Input vs time arrays
Telec_f = NaN(n_t,1);
A_vt_f = NaN(n_t,1);

Q_hyd_f = NaN(n_t,1);

for t_f = 1:1:n_t-1

    J_f = reshape(J_cost(t_f+1,:), n_P1, n_P2, n_theta_dot);

    % All the forward dynamic calculations
    Q_f = Cd*A_vt_options*sqrt((2/rho)*abs(P_H-P2_f(t_f)))*sign(P_H-P2_f(t_f));

    Q_hyd_f(t_f) = (hyd_D/(2*pi))*theta_dot_f(t_f);

    P2dot_f = (beta/(V2_0))*(Q_f - Q_hyd_f(t_f));
    P2next_f = P2_f(t_f) + P2dot_f * dt;

    P1dot_f = (beta/(V1_0)) * (Q_hyd_f(t_f));
    P1next_f = (P1_f(t_f) + P1dot_f * dt) * ones(1,U_input);

    hyd_torque_f = (hyd_D/(2*pi))*(P2_f(t_f) - P1_f(t_f));
    theta_doubledot_f = (1/(J_hyd + J_elec))*(hyd_torque_f - Telec_options);
    theta_dot_next_f = theta_dot_f(t_f) + theta_doubledot_f * dt;
    theta_next_f = ((P1next_f-P_M)*V1_0*2*pi)/(beta*hyd_D);

    cost_next_f = interpn(P1, P2, theta_dot, J_f, P1next_f,...    
    P2next_f,theta_dot_next_f,'linear');

    cost_next_f(isnan(cost_next_f)) = max(J_f(:));

   % Cost and energy calculations

   % Throttling losses
   throttling_now_forward = Q_f * (P_H-P2_f(t_f));

   Qnext_f = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-P2next_f))...
        .*sign(P_H-P2next_f);

   throttling_next_f = Qnext_f .* (P_H-P2next_f);

   energy_throttle_loss_f = 0.5*(throttling_now_forward +...
        throttling_next_f)*dt;

   total_cost_f = cost_next_f + energy_throttle_loss_f;

   % Cost minimization
   % The heart of the forward path planning
   [J_cost_f(t_f),u_opt_f(t_f)] = min(total_cost_f, [], 'omitnan');

   % Setting up the optimal path
   P2_f(t_f+1) = P2next_f(u_opt_f(t_f));
   P1_f(t_f+1) = P1next_f(u_opt_f(t_f));

   Q_forward(t_f) = Q_f(u_opt_f(t_f));
   Q_next_forward(t_f) = Qnext_f(u_opt_f(t_f));

   theta_dot_f(t_f+1) = theta_dot_next_f(u_opt_f(t_f));
   theta_f(t_f+1) = theta_next_f(u_opt_f(t_f));

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
plot(t, theta_f, 'LineWidth',2);
xlabel('Time (s)')
ylabel('m')
grid
title('Displacement')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

sgtitle('Euler Method')

%% Energy plots
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
work_out = cumsum(P1_f.*Q_hyd_f*dt);
plot(t, work_out, 'LineWidth', 3);
xlabel('Time (s)')
ylabel('Energy (J)')
grid
title('Work Out')

subplot(3,2,5)
compress_energy = cumsum(P2_f.*(Q_forward-Q_hyd_f)*dt);
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

sgtitle('Euler Method')

savefig(plot_fig,'Rotary_dp_euler_normal')

