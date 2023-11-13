%%清理工作区数据
clear;                %清除所有变量
close all;                %清图
clc;
% 定义常数 a 和 k
a=100;
k=2;
c=0.2:0.1:0.6;
cmax=max(c);
RFMAX=zeros(1,length(c));
for i=1:1:length(c)
cp=c(i);
rdp=0;
RF=0:0.01:1;
W=zeros(1,length(RF));
CP=zeros(1,length(RF));
for j=1:1:length(RF)
    rdp=RF(j);
% 定义函数
equation = @(w1) (1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rdp)/w1,a,k));
% 提供初始猜测值
initial_guess =0.81; % 可以提供一个合适的初始猜测值

% 使用 fsolve 求解方程
W(j)= fsolve(equation, initial_guess);
if (j~=1)&&(W(j)>W(j-1))
 initial_guess = 0.1; % 可以提供一个合适的初始猜测值
 W(j)= fsolve(equation, initial_guess);
end
CP(j)=cp*(1+RF(j));
if (CP(j)>W(j))
    RFMAX(i)=RF(j-1);
    break;
end
end
Scheme(i).WW_sc=W;
Scheme(i).CP_sc=CP;
end
%HDP=(1/cp-1)*ones(1,length(RF));
figure(2)
RF=0:0.01:RFMAX(1);
plot(RF,Scheme(1).WW_sc(1:length(RF)),'r-','DisplayName',sprintf('c=%0.1f',c(1)),'LineWidth',2);
hold on
RF=0:0.01:RFMAX(1);
plot(RF,Scheme(1).CP_sc(1:length(RF)),'r-','DisplayName',sprintf('w1=%0.1f',c(1)),'LineWidth',2);
for i=2:1:length(c)
color=[0.1*i,1-0.1*i,1-0.1*i];
RF=0:0.01:RFMAX(i);
plot(RF, Scheme(i).WW_sc(1:length(RF)),'color', color,'DisplayName',sprintf('c=%0.1f',c(i)),'LineWidth',2);
RF=0:0.01:RFMAX(i);
plot(RF,Scheme(i).CP_sc(1:length(RF)),'color', color,'DisplayName',sprintf('w1=%0.1f',c(i)),'LineWidth',2);
end
%plot(HDP,W,'bo-','DisplayName','rf=1/cp-1');
xlim([0, 1]); 
ylim([0, 1]);
xlabel('rf');
ylabel('w1');
legend
  




