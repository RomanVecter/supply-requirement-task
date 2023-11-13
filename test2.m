%%清理工作区数据
clear;                %清除所有变量
close all;                %清图
clc;
% 定义常数 a 和 k
a=100;
k=2;
c=0.2:0.1:0.6;
cmax=max(c);
for j=1:1:length(c)
cp=c(j);
equation = @(w0) (1/w0-1)*cp/w0/h_PT(F_fan_inverse(cp/w0,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w0,a,k));
% 提供初始猜测值
initial_guess = 0.1; % 可以提供一个合适的初始猜测值
% 使用 fsolve 求解方程
w0= fsolve(equation, initial_guess);
SchemeW(j).WW_sc=w0:0.01:1;
W=SchemeW(j).WW_sc;
RF=zeros(1,length(W));
CP=zeros(1,length(W));
for i=1:1:length(W)
    w1=W(i);
% 定义函数
equation = @(rf0) (1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k));
% 提供初始猜测值
initial_guess = 0.1; % 可以提供一个合适的初始猜测值
% 使用 fsolve 求解方程
RF(i)= fsolve(equation, initial_guess);
CP(i)=cp*(1+RF(i));
end
Scheme(j).RF_sc=RF;
Scheme(j).CP_sc=CP;
end
W=0:0.01:1;
plot(Scheme(1).RF_sc,SchemeW(1).WW_sc,'r-','DisplayName',sprintf('c=%0.1f',c(1)),'LineWidth',2);

xlim([0, 1]); 
ylim([0,1]); 
hold on
plot(Scheme(1).RF_sc,Scheme(1).CP_sc,'r-','DisplayName',sprintf('w1=%0.1f',c(1)),'LineWidth',2);
for i=2:1:length(c)
color=[0.1*i,1-0.1*i,1-0.1*i];
plot(Scheme(i).RF_sc,SchemeW(i).WW_sc,'color', color,'DisplayName',sprintf('c=%0.1f',c(i)),'LineWidth',2);
plot(Scheme(i).RF_sc,Scheme(i).CP_sc,'r-','color', color,'DisplayName',sprintf('w1=%0.1f',c(i)),'LineWidth',2);
end
% plot(HDP,W,'bo-','DisplayName','rf=1/cmax-1');

xlabel('rf');
ylabel('W');
legend