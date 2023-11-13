%%model2
%%清理工作区数据
clear;                %清除所有变量
close all;                %清图
clc;
%%参数定义
%%供应商
%%制造商
%%零售商
%%韦伯分布参数控制
a=100;%%缩放因子
k=2;%%形状参数 
%%隐函数方程组求解
c=0.2:0.05:0.8;
rf0_solution=zeros(1,length(c));
w1_solution=zeros(1,length(c));
for i=1:1:length(c)
    cp=c(i);
F=@(rf0,w1) cp*(1+rf0)*(F_fan_inverse(cp,a,k)-F_fan_inverse(cp*(1+rf0)/w1,a,k))+(1-w1)*integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+w1*integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k))-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp,a,k));
G=@(rf0,w1) (1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k));
% 初始猜测值
initial_guess = [1.0, 1.0];
equations = @(x) [F(x(1), x(2)); G(x(1), x(2))];
%options = optimoptions('fsolve', 'Display', 'iter');  % 设置显示迭代过程
solution = fsolve(equations,initial_guess);
rf0_solution(i)= solution(1);
w1_solution(i) = solution(2);
end
figure(1)
xlim([min(c), max(c)]); 
plot(c,rf0_solution,'r*-','DisplayName','c-rf0曲线');
hold on
xlim([min(c), max(c)]); 
plot(c,w1_solution,'bo-','DisplayName','c-w1曲线');
legend
clc;
%代数求解
%计算数据
c=0.2:0.1:0.6;
for i=1:1:length(c)
    cp=c(i);
    [cp,rf0,w1,RFMAX,WT,QT,KT,ST,MT,LT]=Model2_FUC(cp,a,k);
    Scheme(i).c_sc=cp;
    Scheme(i).rf0_sc=rf0;
    Scheme(i).w1_sc=w1;
    Scheme(i).RFMAX_sc=RFMAX;
    Scheme(i).w_sc=WT;
    Scheme(i).Q_sc=QT;
    Scheme(i).K_sc=KT;
    Scheme(i).S_sc=ST;
    Scheme(i).M_sc=MT;
    Scheme(i).L_sc=LT;
end
%%作图
figure(2)
set(gcf, 'Position', [480, 100, 600, 600]);%%图像居中
subplot(3,1,1)
rf=0.0:0.01:Scheme(1).RFMAX_sc;
plot(rf, Scheme(1).w_sc,'r-','DisplayName',sprintf('c=%0.1f,RFMAX=%0.2f,w1=%0.2f',Scheme(1).c_sc,Scheme(1).RFMAX_sc,Scheme(1).w1_sc),'LineWidth',2);
hold on
for i=2:1:length(c)
color=[0.1*i,1-0.1*i,1-0.1*i];
rf=0.0:0.01:Scheme(i).RFMAX_sc;
plot(rf, Scheme(i).w_sc,'color', color,'DisplayName',sprintf('c=%0.1f,RFMAX=%0.2f,w1=%0.2f',Scheme(i).c_sc,Scheme(i).RFMAX_sc,Scheme(i).w1_sc),'LineWidth',2);
end
xlim([0, 1]); 
% ylim([0.6, 1.2]);
xlabel('rf');
ylabel('w');
title('a)rf-w曲线');
legend
subplot(3,1,2)
rf=0.0:0.01:Scheme(1).RFMAX_sc;
plot(rf, Scheme(1).Q_sc,'r-','DisplayName',sprintf('c=%0.1f,RFMAX=%0.2f,w1=%0.2f',Scheme(1).c_sc,Scheme(1).RFMAX_sc,Scheme(1).w1_sc),'LineWidth',2);
hold on
for i=2:1:length(c)
color=[0.1*i,1-0.1*i,1-0.1*i];
rf=0.0:0.01:Scheme(i).RFMAX_sc;
plot(rf, Scheme(i).Q_sc,'color', color,'DisplayName',sprintf('c=%0.1f,RFMAX=%0.2f,w1=%0.2f',Scheme(i).c_sc,Scheme(i).RFMAX_sc,Scheme(i).w1_sc),'LineWidth',2);
end
xlim([0, 1]); 
% ylim([0, 200]);
xlabel('rf');
ylabel('Q');
title('b)rf-Q曲线');
legend
subplot(3,1,3)
rf=0.0:0.01:Scheme(1).RFMAX_sc;
plot(rf, Scheme(1).K_sc,'r-','DisplayName',sprintf('c=%0.1f,RFMAX=%0.2f,w1=%0.2f',Scheme(1).c_sc,Scheme(1).RFMAX_sc,Scheme(1).w1_sc),'LineWidth',2);
hold on
for i=2:1:length(c)
color=[0.1*i,1-0.1*i,1-0.1*i];
rf=0.0:0.01:Scheme(i).RFMAX_sc;
plot(rf, Scheme(i).K_sc,'color', color,'DisplayName',sprintf('c=%0.1f,RFMAX=%0.2f,w1=%0.2f',Scheme(i).c_sc,Scheme(i).RFMAX_sc,Scheme(i).w1_sc),'LineWidth',2);
end
xlim([0, 1]); 
% ylim([0, 50]);
xlabel('rf');
ylabel('K');
title('c)rf-K曲线');
legend
%%均衡利润曲线
figure(3)
set(gcf, 'Position', [480, 100, 600, 600]);%%图像居中
subplot(3,1,1)
rf=0.0:0.01:Scheme(1).RFMAX_sc;
plot(rf, Scheme(1).S_sc,'r-','DisplayName',sprintf('c=%0.1f,RFMAX=%0.2f,w1=%0.2f',Scheme(1).c_sc,Scheme(1).RFMAX_sc,Scheme(1).w1_sc),'LineWidth',2);
hold on
for i=2:1:length(c)
color=[0.1*i,1-0.1*i,1-0.1*i];
rf=0.0:0.01:Scheme(i).RFMAX_sc;
plot(rf, Scheme(i).S_sc,'color', color,'DisplayName',sprintf('c=%0.1f,RFMAX=%0.2f,w1=%0.2f',Scheme(i).c_sc,Scheme(i).RFMAX_sc,Scheme(i).w1_sc),'LineWidth',2);
end
xlim([0, 1]); 
% ylim([0, 60]);
xlabel('rf');
ylabel('S');
title('a)rf-S曲线');
legend
subplot(3,1,2)
rf=0.0:0.01:Scheme(1).RFMAX_sc;
plot(rf, Scheme(1).M_sc,'r-','DisplayName',sprintf('c=%0.1f,RFMAX=%0.2f,w1=%0.2f',Scheme(1).c_sc,Scheme(1).RFMAX_sc,Scheme(1).w1_sc),'LineWidth',2);
hold on
for i=2:1:length(c)
color=[0.1*i,1-0.1*i,1-0.1*i];
rf=0.0:0.01:Scheme(i).RFMAX_sc;
plot(rf, Scheme(i).M_sc,'color', color,'DisplayName',sprintf('c=%0.1f,RFMAX=%0.2f,w1=%0.2f',Scheme(i).c_sc,Scheme(i).RFMAX_sc,Scheme(i).w1_sc),'LineWidth',2);
end
xlim([0, 1]); 
% ylim([0, 100]);
xlabel('rf');
ylabel('M');
title('b)rf-M曲线');
legend
subplot(3,1,3)
rf=0.0:0.01:Scheme(1).RFMAX_sc;
plot(rf, Scheme(1).L_sc,'r-','DisplayName',sprintf('c=%0.1f,RFMAX=%0.2f,w1=%0.2f',Scheme(1).c_sc,Scheme(1).RFMAX_sc,Scheme(1).w1_sc),'LineWidth',2);
hold on
for i=2:1:length(c)
    color=[0.1*i,1-0.1*i,1-0.1*i];
    rf=0.0:0.01:Scheme(i).RFMAX_sc;
plot(rf, Scheme(i).L_sc,'color', color,'DisplayName',sprintf('c=%0.1f,RFMAX=%0.2f,w1=%0.2f',Scheme(i).c_sc,Scheme(i).RFMAX_sc,Scheme(i).w1_sc),'LineWidth',2);
end
xlim([0, 1]); 
% ylim([40, 100]);
xlabel('rf');
ylabel('L');
title('c)rf-L曲线');
legend


