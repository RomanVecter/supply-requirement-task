%%清理工作区数据
clear;                %清除所有变量
close all;                %清图
clc;
% 定义常数 a 和 k
a=100;
k=2;
RFMAX=0;
cp=0.4;
rdp=0;
RF=0:0.01:1;
W=zeros(1,length(RF));
CP=zeros(1,length(RF));
for j=1:1:length(RF)
    rdp=RF(j);
% 定义函数
equation = @(w1) (1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rdp)/w1,a,k));
% 提供初始猜测值
initial_guess =0.82; % 可以提供一个合适的初始猜测值

% 使用 fsolve 求解方程
W(j)= fsolve(equation, initial_guess);
if (j~=1)&&(W(j)>W(j-1))
 initial_guess = 0.1; % 可以提供一个合适的初始猜测值
 W(j)= fsolve(equation, initial_guess);
end
CP(j)=cp*(1+RF(j));
if (CP(j)>W(j))
    RFMAX=RF(j-1);
    break;
end
end
RF=0:0.01:RFMAX;
plot(RF,W(1:length(RF)),'r-','DisplayName',sprintf('c=%0.1f',cp),'LineWidth',2);
hold on
RF=0:0.01:RFMAX(1);
plot(RF,CP(1:length(RF)),'r-','DisplayName',sprintf('w1=%0.1f',cp),'LineWidth',2);
% xlim([0,RFMAX]);
% ylim([0,1]);