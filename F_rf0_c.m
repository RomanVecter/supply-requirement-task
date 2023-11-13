function [rf0]=F_rf0_c(c,a,k)
equation = @(rf0) F_fan((F_fan(H_inverse(rf0/(1+rf0),a,k),a,k)/c*H_inverse(rf0/(1+rf0),a,k)),a,k)/c-1-rf0;
% 提供初始猜测值
initial_guess = 0.5; % 可以提供一个合适的初始猜测值
% 使用 fsolve 求解方程
rf0= fsolve(equation, initial_guess);