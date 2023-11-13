%%model5
%%验证供应链利润与制造商利润曲线图
%%修改15行的c，即可实现不同c对应的利润出图，
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
%%自己设定c
c=0.3;

[cp,rf0,w2,w1_min,RFMAX,L_F,M_F,R_F,L_B,M_B,R_B]=Model5_FUC(c,a,k);
L_limit=0;
H_limit=100;

%%
x1=rf0;
x2=RFMAX;
x_limit=x2;
figure(1)
set(gcf, 'Position', [480, 100, 1000, 600]);%%图像居中
rf=0.01:0.01:1;
plot(x1* ones(size(rf)), rf*100,'r--','DisplayName',sprintf('rf=rf0,    c=%0.1f,    rf0=%0.2f',c,rf0),'LineWidth', 2 );  % 绘制红色直线
hold on
plot(x2* ones(size(rf)), rf*100, 'b--','DisplayName',sprintf('rf=RFMAX,     c=%0.1f,    RFMAX=%0.2f',c,x2),'LineWidth', 2);  % 绘制蓝色直线
xlim([0, x_limit*1.2]); 
rf=0.01:0.01:RFMAX;
plot(rf, L_F,'color', [0.1,0.9,0.5],'DisplayName',sprintf('供应链利润曲线L_F,c=%0.1f*(1+RFMAX=%0.2f)<w1=%0.2f',c,RFMAX,w1_min),'LineWidth',2); 
plot(rf, M_F,'color', [0.5,0.5,0.3],'DisplayName',sprintf('供应链利润曲线M_F,c=%0.1f*(1+RFMAX=%0.2f)<w1=%0.2f',c,RFMAX,w1_min),'LineWidth',2);
plot(rf, R_F,'color', [0.8,0.6,0.2],'DisplayName',sprintf('供应链利润曲线R_F,c=%0.1f*(1+RFMAX=%0.2f)<w1=%0.2f',c,RFMAX,w1_min),'LineWidth',2); 

xlabel('rf');
ylabel('供应链利润');
title('供应链利润曲线');
legend

figure(2)
set(gcf, 'Position', [480, 100, 1000, 600]);%%图像居中
rf=0.01:0.01:1;
plot(x1* ones(size(rf)), rf*100,'r--','DisplayName',sprintf('rf=rf0,    c=%0.1f,    rf0=%0.2f',c,rf0),'LineWidth', 2 );  % 绘制红色直线
hold on
plot(x2* ones(size(rf)), rf*100, 'b--','DisplayName',sprintf('rf=RFMAX,     c=%0.1f,    RFMAX=%0.2f',c,x2),'LineWidth', 2);  % 绘制蓝色直线
xlim([0, x_limit*1.2]); 
rf=0.01:0.01:RFMAX;
plot(rf, L_B,'color', [0.1,0.9,0.5],'DisplayName',sprintf('制造商利润曲线L_B,c=%0.1f*(1+RFMAX=%0.2f)<w1=%0.2f',c,RFMAX,w1_min),'LineWidth',2);  
plot(rf, M_B,'color', [0.5,0.5,0.3],'DisplayName',sprintf('制造商利润曲线M_B,c=%0.1f*(1+RFMAX=%0.2f)<w1=%0.2f',c,RFMAX,w1_min),'LineWidth',2);
plot(rf, R_B,'color', [0.8,0.6,0.2],'DisplayName',sprintf('制造商利润曲线R_B,c=%0.1f*(1+RFMAX=%0.2f)<w1=%0.2f',c,RFMAX,w1_min),'LineWidth',2); 

xlabel('rf');
ylabel('制造商利润');
title('制造商利润曲线');
legend