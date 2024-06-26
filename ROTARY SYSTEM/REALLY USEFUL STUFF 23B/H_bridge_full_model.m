% Clearing everything

clc; clear; close all;

% Defining the symbols. Note that P1, Q

syms Q1 Q2 Q3 Q4 Q5 k_A k_B k_C k_D k Pdes P1 P3 P4 alpha_1 alpha_2

% Setting up the simultaneous equations

eqn(1) = Q1 == k_A*sqrt(Pdes-P3);
eqn(2) = Q2 == k_B*sqrt(P3-P1);
eqn(3) = Q2 == Q1 - Q3;
eqn(4) = Q4 == k_C*sqrt(P4-P1);
eqn(5) = Q5 == k_D*sqrt(Pdes-P4);
eqn(6) = Q4 == Q3 + Q5;

% Solving the simultaneous equations

S = solve(eqn, [Q1, Q2, Q4, Q5, P3, P4]);

disp(S.Q1)
disp(S.Q2)
disp(S.Q4)
disp(S.Q5)
disp(S.P3)
disp(S.P4)