%%清理工作区数据
clear;                %清除所有变量
close all;                %清图
clc;
% 定义常数 a 和 k
a=100;
k=2;
cp=0.5;
rf0=0;
equation = @(w1) (1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k));
% 提供初始猜测值
initial_guess =1; % 可以提供一个合适的初始猜测值
% 使用 fsolve 求解方程

w1= fsolve(equation, initial_guess);
if isreal(w1)
    disp('w1这是一个实数。');
else
    disp('w1这是一个复数。');
end
cp*(1+rf0)/w1

F_fan_inverse(cp*(1+rf0)/w1,a,k)
(1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)
integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))
integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k))
FP=(1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k))