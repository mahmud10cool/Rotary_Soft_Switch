% Author: Mahmud Suhaimi Ibrahim
% Dynamic Programming Code

%% Dynamic Programming for Soft Switch
clc; clear;
format long;

% State, Input and Time array sizes
P1_size = 11;
P2_size = 11;
Xdot_size = 251;
Av_size = 1;
Felec_size = 1;
FINALTIME = 60e-3;
STEPSIZE = 1e-5;
TIMESIZE = 1 + FINALTIME/STEPSIZE;
weight = 5;

save ARRAY.mat P1_size P2_size Xdot_size Av_size Felec_size FINALTIME
arr_size = [num2str(P1_size),' ',num2str(P2_size),' ', ...
    num2str(Xdot_size),' ', num2str(Av_size),' ',num2str(Felec_size),...
    ' ', num2str(FINALTIME)];
disp(arr_size)

% Added Mass
delm = 0;

% Maximum and minimum areas
Ap = 0.25*pi*(20e-3)^2;              % Plunger cross-sectional area
Av_max = 0.5*Ap;

% Maximum and minimum forces
Fmin = 0;
Fmax = 0;

% Maximum and minimum velocities
Xdot_min = -500;
Xdot_max = 500;


%% Parameter values
% Plunger
Dp = 5e-3;                     % Plunger diameter
Hp = 20e-3;                     % Plunger height
Ap = 0.25*pi*Dp^2;              % Plunger cross-sectional area
rho_ss = 7800;                  % Density of steel
mass = rho_ss*Ap*Hp;            % Mass of the plunger
V2_0 = 20e-6;                   % Initial volume below plunger                  
c = 0;

% Hydraulic oil
beta = 1.8e9;                   % Bulk Modulus
rho = 870;                      % Density of fluid

% Linear actuator chamber
V1_0 = 2e-3;                    % Initial volume in linear actuator

% Valve properties
Cd = 0.6;                       % Coefficient of discharge

P_M = 10e6;                     % Medium pressure rail
P_H = 20e6;                     % High pressure rail

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
n_Xdot = Xdot_size;
Xdot = linspace(Xdot_min, Xdot_max, n_Xdot);

if ~any(Xdot==0)
    return
end

[P1_matrix, P2_matrix, Xdot_matrix] = ndgrid(P1,P2,Xdot);

P1_states = P1_matrix(:);
P2_states = P2_matrix(:);
Xdot_states = Xdot_matrix(:);

% Number of discretizations by changing ndgrid matrix into column vector
States = length(P1_states);

%% Input discretization

n_Av = Av_size;                       
A_vt = linspace(0,Av_max,n_Av); 

n_Felec = Felec_size;
Felec = linspace(Fmin,Fmax,n_Felec); 

[A_vt_matrix,Felec_matrix] = ndgrid(A_vt,Felec);

A_vt_options = A_vt_matrix(:)';
Felec_options = Felec_matrix(:)';

U_input = length(A_vt_options); 

%% Dynamic programming

% Initialising the cost matrix
J_cost = NaN(n_t,States);

% This is the thing altering the results the most
ind_final_P1 = find(P1 == P_H);
ind_final_P2 = find(P2 == P_H);
ind_final_Xdot = find(Xdot == 0);

ind_final = sub2ind(size(P1_matrix), ind_final_P1, ind_final_P2, ind_final_Xdot);

%  Cost at final time
J_cost(end,:) = weight*(1e-5*abs(P1_states-P1(ind_final_P1)).^1 + 1e-5*abs(P2_states - P2(ind_final_P2)).^1 + ...
    10*abs((Xdot_states-Xdot(ind_final_Xdot))).^1);
% J_cost(end,:) = 0;

% Initialising the optimal input (index) matrix
u_opt_matrix = NaN(n_t,States); 
u_opt_matrix(end,:) = 1;

%% Precalculations

% I will use 'k' for P2, 'l' for P1 and 'm' for Xdot

%% First term
Q1 = Cd.*A_vt_options.*sqrt((2/rho)*abs(P_H-P2_states)).*sign(P_H-P2_states);
    
X1 = (V1_0/Ap)*(1-exp((P_M-P1_states)/beta));
    
P2dot = (beta./(V2_0+Ap*X1)).*(Q1-Ap*Xdot_states);

P1dot = (beta./(V1_0-Ap*X1)).*(Ap*Xdot_states);

Xdoubledot = (1/mass) * ((P2_states - P1_states)*Ap - Felec_options);

k1 = P2dot;
l1 = P1dot;
m1 = Xdoubledot;

%% Second term
Q2 = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-(P2_states+0.5*dt*k1))).*sign(P_H-(P2_states+0.5*dt*k1));

X2 = (V1_0/Ap)*(1-exp((P_M-(P1_states+0.5*dt*l1))/beta));

k2 = (beta./(V2_0+Ap*X2)).*(Q2-Ap*(Xdot_states+0.5*dt*m1));
l2 = (beta./(V1_0-Ap*X2)).*(Ap*(Xdot_states+0.5*dt*m1));
m2 = (1/mass) * (((P2_states+0.5*dt*k1) - (P1_states+0.5*dt*l1))*Ap - Felec_options);

%% Third term
Q3 = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-(P2_states+0.5*dt*k2))).*sign(P_H-(P2_states+0.5*dt*k2));

X3 = (V1_0/Ap)*(1-exp((P_M-(P1_states+0.5*dt*l2))/beta));

k3 = (beta./(V2_0+Ap*X3)).*(Q3-Ap*(Xdot_states+0.5*dt*m2));
l3 = (beta./(V1_0-Ap*X3)).*(Ap*(Xdot_states+0.5*dt*m2));
m3 = (1/mass) * (((P2_states+0.5*dt*k2) - (P1_states+0.5*dt*l2))*Ap - Felec_options);

%% Fourth term
Q4 = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-(P2_states+dt*k3))).*sign(P_H-(P2_states+dt*k3));

X4 = (V1_0/Ap)*(1-exp((P_M-(P1_states+dt*l3))/beta));

k4 = (beta./(V2_0+Ap*X4)).*(Q4-Ap*(Xdot_states+dt*m3));
l4 = (beta./(V1_0-Ap*X4)).*(Ap*(Xdot_states+dt*m3));
m4 = (1/mass) * (((P2_states+dt*k3) - (P1_states+dt*l3))*Ap - Felec_options);

%% Integration
P2next = P2_states + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
P1next = P1_states + (dt/6)*(l1 + 2*l2 + 2*l3 + l4);
Xdotnext = Xdot_states + (dt/6)*(m1 + 2*m2 + 2*m3 + m4);

%% Throttling loss calculations
throttling_1 = Q1 .* (P_H-P2_states);

throttling_2 = Q2 .* (P_H-(P2_states+0.5*dt*k1));

throttling_3 = Q3 .* (P_H-(P2_states+0.5*dt*k2));

throttling_4 = Q4 .* (P_H-(P2_states+dt*k3));

energy_throttle_loss = 0.5*(throttling_1+2*throttling_2+2*throttling_3+throttling_4)*dt;

%% Starting the DP

startTimer = tic;

disp('DP started successfully.')

for t_ind = n_t-1:-1:1

    J = reshape(J_cost(t_ind+1,:), n_P1, n_P2, n_Xdot);

    cost_next = interpn(P1,P2,Xdot,J,P1next,...
            P2next,Xdotnext,'linear'); 
        
    cost_next(isnan(cost_next)) = max(J(:));

    total_loss = cost_next + energy_throttle_loss; 

   [J_cost(t_ind,:),u_opt_matrix(t_ind,:)] = min(total_loss, ...
         [], 2,'omitnan'); 
    
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
P2_f = NaN(n_t, 1);
P1_f = NaN(n_t, 1);
Xdot_f = NaN(n_t, 1);

% Flow arrays
Q_forward = NaN(n_t,1);
Q_next_forward = NaN(n_t,1);

% Find the right initial states
ind_P1 = find(P1 == P_M);
ind_P2 = find(P2 == P_M);
ind_Xdot = find(Xdot == 0);

% Initial states
P1_f(1,1) = P1(ind_P1);
P2_f(1,1) = P2(ind_P2);
Xdot_f(1,1) = Xdot(ind_Xdot);

X1_f = NaN(n_t, 1);
X2_f = NaN(n_t, 1);
X3_f = NaN(n_t, 1);
X4_f = NaN(n_t, 1);

J_cost_f = NaN(n_t, 1);
u_opt_f = NaN(n_t, 1);
loss_f = NaN(n_t,1);
total_loss_f = NaN(n_t,1);
total_loss_f(1) = 0;

% Input vs time arrays
Felec_f = NaN(n_t,1);
A_vt_f = NaN(n_t,1);

for t_f = 1:n_t-1 
 
    J_f = reshape(J_cost(t_f+1,:), n_P1, n_P2, n_Xdot);

    %% First term
    
    Q1_f = Cd*A_vt_options*sqrt((2/rho)*abs(P_H-P2_f(t_f)))*sign(P_H-P2_f(t_f));
    
    X1_f(t_f) = (V1_0/Ap)*(1-exp((P_M-P1_f(t_f))/beta));
    
    P2dot_f = (beta/(V2_0+Ap*X1_f(t_f)))*(Q1_f-Ap*Xdot_f(t_f));

    P1dot_f = (beta/(V1_0-Ap*X1_f(t_f)))*(Ap*Xdot_f(t_f));

    Xdoubledot_f = (1/mass) * ((P2_f(t_f) - P1_f(t_f))*Ap - Felec_options);

    k1_f = P2dot_f;
    l1_f = P1dot_f;
    m1_f = Xdoubledot_f;

    %% Second term
    Q2_f = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-(P2_f(t_f)+0.5*dt*k1_f))).*sign(P_H-(P2_f(t_f)+0.5*dt*k1_f));

    X2_f = (V1_0/Ap)*(1-exp((P_M-(P1_f(t_f)+0.5*dt*l1_f))/beta));

    k2_f = (beta./(V2_0+Ap*X2_f)).*(Q2_f-Ap*(Xdot_f(t_f)+0.5*dt*m1_f));
    l2_f = (beta./(V1_0-Ap*X2_f)).*(Ap*(Xdot_f(t_f)+0.5*dt*m1_f));
    m2_f = (1/mass) * (((P2_f(t_f)+0.5*dt*k1_f) - (P1_f(t_f)+0.5*dt*l1_f))*Ap - Felec_options);

    %% Third term
    Q3_f = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-(P2_f(t_f)+0.5*dt*k2_f))).*sign(P_H-(P2_f(t_f)+0.5*dt*k2_f));

    X3_f = (V1_0/Ap)*(1-exp((P_M-(P1_f(t_f)+0.5*dt*l2_f))/beta));

    k3_f = (beta./(V2_0+Ap*X3_f)).*(Q3_f-Ap*(Xdot_f(t_f)+0.5*dt*m2_f));
    l3_f = (beta./(V1_0-Ap*X3_f)).*(Ap*(Xdot_f(t_f)+0.5*dt*m2_f));
    m3_f = (1/mass) * (((P2_f(t_f)+0.5*dt*k2_f) - (P1_f(t_f)+0.5*dt*l2_f))*Ap - Felec_options);

    %% Fourth term
    Q4_f = Cd*A_vt_options.*sqrt((2/rho)*abs(P_H-(P2_f(t_f)+dt*k3_f))).*sign(P_H-(P2_f(t_f)+dt*k3_f));

    X4_f = (V1_0/Ap)*(1-exp((P_M-(P1_f(t_f)+dt*l3_f))/beta));

    k4_f = (beta./(V2_0+Ap*X4_f)).*(Q4_f-Ap*(Xdot_f(t_f)+dt*m3_f));
    l4_f = (beta./(V1_0-Ap*X4_f)).*(Ap*(Xdot_f(t_f)+dt*m3_f));
    m4_f = (1/mass) * (((P2_f(t_f)+dt*k3_f) - (P1_f(t_f)+dt*l3_f))*Ap - Felec_options);

    %% Integration
    P2next_f = P2_f(t_f) + (dt/6)*(k1_f + 2*k2_f + 2*k3_f + k4_f);
    P1next_f = P1_f(t_f) + (dt/6)*(l1_f + 2*l2_f + 2*l3_f + l4_f);
    Xdotnext_f = Xdot_f(t_f) + (dt/6)*(m1_f + 2*m2_f + 2*m3_f + m4_f);

    %% X
    Xnext_f = (V1_0/Ap)*(1-exp((P_M-P1next_f)/beta));

    %% Cost analysis
    cost_next_f = interpn(P1, P2, Xdot, J_f, P1next_f,...
        P2next_f,Xdotnext_f,'linear');

    cost_next_f(isnan(cost_next_f)) = max(J(:));

    %% Throttling losses

    throttling_1_f = Q1 .* (P_H-P2_f(t_f));

    throttling_2_f = Q2 .* (P_H-(P2_f(t_f)+0.5*dt*k1_f));

    throttling_3_f = Q3 .* (P_H-(P2_f(t_f)+0.5*dt*k2_f));

    throttling_4_f = Q4 .* (P_H-(P2_f(t_f)+dt*k3_f));

    energy_throttle_loss_f = 0.5*(throttling_1_f+throttling_2_f+throttling_3_f+throttling_4)*0.5*dt;

    total_cost_f = cost_next_f + energy_throttle_loss_f;
 
    % Cost minimization
    % The heart of the forward path planning
    [J_cost_f(t_f),u_opt_f(t_f)] = min(total_cost_f, [], 'omitnan');
    
    % Setting up the optimal path
    P2_f(t_f+1) = P2next_f(u_opt_f(t_f));
    P1_f(t_f+1) = P1next_f(u_opt_f(t_f));
    
    Q_forward(t_f) = Q1_f(u_opt_f(t_f));
    Q_next_forward(t_f) = Qnext_f(u_opt_f(t_f));
    
    Xdot_f(t_f+1) = Xdotnext_f(u_opt_f(t_f));
    X1_f(t_f+1) = Xnext_f(u_opt_f(t_f));
    
    A_vt_f(t_f) = A_vt_options(u_opt_f(t_f));
    Felec_f(t_f) = Felec_options(u_opt_f(t_f));
    
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
plot(t, Felec_f, 'LineWidth',2);
xlabel('Time (s)')
ylabel('Force (N)')
grid
title('Electric Force Applied')

subplot(3,2,4)
plot(t, A_vt_f, 'LineWidth',2);
xlabel('Time (s)')
ylabel('Area (m^2)')
grid
title('Valve Area Opening')
ylim([0 1.5*max(A_vt_f)])

subplot(3,2,5)
plot(t, Xdot_f, 'LineWidth',2);
xlabel('Time (s)')
ylabel('ms^-1')
grid
title('Velocity')

subplot(3,2,6)
plot(t, X1_f, 'LineWidth',2);
xlabel('Time (s)')
ylabel('m')
grid
title('Displacement')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

% Flow Plots
plot_fig(2) = figure;

subplot(2,1,1)
plot(t, Q_forward, Linewidth=2)
xlabel('Time (s)')
ylabel('Flow (m^3s^-1)')
grid
title('Flow Rate')

subplot(2,1,2)
plot(t, Q_next_forward, Linewidth=2)
xlabel('Time (s)')
ylabel('Flow (m^3s^-1)')
grid
title('Flow Rate')

% Energy Plots
plot_fig(3) = figure;

subplot(3,2,1)
plot(t, total_loss_f, 'LineWidth', 3);
xlabel('Time (s)')
ylabel('Energy (J)')
grid
title('Energy Loss')

subplot(3,2,2)
regen = cumsum(Felec_f.*Xdot_f*dt);
plot(t, regen, 'LineWidth', 3);
xlabel('Time (s)')
ylabel('Energy (J)')
grid
title('Regenerated')

subplot(3,2,3)
plot(t, cumsum(P_H.*Q_forward*dt), 'LineWidth', 3);
xlabel('Time (s)')
ylabel('Energy (J)')
grid
title('Energy In')

subplot(3,2,4)
work_out = cumsum(P1_f.*Xdot_f*Ap*dt);
plot(t, work_out, 'LineWidth', 3);
xlabel('Time (s)')
ylabel('Energy (J)')
grid
title('Work Out')

subplot(3,2,5)
compress_energy = cumsum(P2_f.*(Q_forward-Ap*Xdot_f)*dt);
plot(t, compress_energy, 'LineWidth', 3);
xlabel('Time (s)')
ylabel('Energy (J)')
grid
title('Compressed Energy')

subplot(3,2,6)
kinetic_energy = 0.5*mass*Xdot_f.^2;
energy_out = total_loss_f + regen + work_out + compress_energy + kinetic_energy;
plot(t, energy_out, 'LineWidth', 3);
xlabel('Time (s)')
ylabel('Energy (J)')
grid
title('Energy Out')

set(findall(gcf,'-property','FontSize'),'FontSize',15)

savefig(plot_fig,'Three_figures.fig')

% Recovered energy percentage
final_total = total_loss_f(end-1) + regen(end-1);
percent_recovered = regen(end-1)/final_total;
dsp_recover = ['Percentage Recovered: ', num2str(percent_recovered*100)];
disp(dsp_recover)

%% Checking something
ind = sub2ind(size(J), ind_P1, ind_P2, ind_Xdot);

back = J_cost(1,ind);

front = total_loss_f(end);

