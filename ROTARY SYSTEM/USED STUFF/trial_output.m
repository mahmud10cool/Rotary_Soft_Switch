clc; clear; close;

tspan = 0:0.01:1;
x0 = [5; 5];

fun = @(t,x) something(t,x);
[t, x] = ode45(@something, tspan, x0);

[~,u] = cellfun(fun,num2cell(t),num2cell(x,2),'uni',0); 
u = cell2mat(u);

figure(1)
plot(t,u)
xlabel('Time (s)')
ylabel('u')
title('u vs t')
grid

figure(2)
plot(t,x(:,1))
xlabel('Time (s)')
ylabel('x1')
title('x1 vs t')
grid

figure(3)
plot(t,x(:,2))
xlabel('Time (s)')
ylabel('x2')
title('x2 vs t')
grid

function [dxdt, u] = something(t,x)

u = t^2;

dxdt = zeros(2,1);

dxdt(1) = x(1);
dxdt(2) = x(1)^2+3 + x(2);

end