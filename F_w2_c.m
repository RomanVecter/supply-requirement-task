function [w2]=F_w2_c(c,a,k)
equation = @(w2) F_fan(w2/c*F_fan_inverse(w2,a,k),a,k)*(1-H(F_fan_inverse(w2,a,k),a,k))/c-1;
% 提供初始猜测值
initial_guess = 0.5; % 可以提供一个合适的初始猜测值
% 使用 fsolve 求解方程
w2= fsolve(equation, initial_guess);