%%model1
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
%%隐函数
syms c rf0
eqn=F_fan((F_fan(H_inverse(rf0/(1+rf0),a,k),a,k)/c*H_inverse(rf0/(1+rf0),a,k)),a,k)/c-1-rf0;
fimplicit(eqn, [0, 1, 0, 1],'DisplayName','rf0'); % 指定 x 和 y 的范围
figure(1);
xlabel('c');
ylabel('rf0');
hold on
c=0:0.01:1;
rf=zeros(1,length(c));
for i=1:1:(length(c))
    rf(i)=1/c(i)-1;
end
plot(c,rf,'DisplayName','rf');
xlim([0, 1]); 
ylim([0, 1]);
legend;
%%c-w2曲线
syms c w2
epn=F_fan(w2/c*F_fan_inverse(w2,a,k),a,k)*(1-H(F_fan_inverse(w2,a,k),a,k))/c-1;
figure(2);
fimplicit(epn, [0, 1, 0, 1],'DisplayName','w2'); % 指定 x 和 y 的范围
hold on
c=0:0.01:1;
w2=c;
plot(c,w2,'DisplayName','w2=c');
xlabel('c');
ylabel('w2');
legend;

%%计算数据
c=0.1:0.1:0.5;
for i=1:1:length(c)
    [cc,rf0,w2,w,Q,q,S,M,R,L]=Model1_FUC(c(i),a,k);
    Scheme(i).c_sc=cc;
    Scheme(i).rf0_sc=rf0;
    Scheme(i).w2_sc=w2;
    Scheme(i).w_sc=w;
    Scheme(i).Q_sc=Q;
    Scheme(i).q_sc=q;
    Scheme(i).S_sc=S;
    Scheme(i).M_sc=M;
    Scheme(i).R_sc=R;
    Scheme(i).L_sc=L;
end
%%作图
rf=0.01:0.01:1;
figure(3)
set(gcf, 'Position', [480, 100, 600, 600]);%%图像居中
subplot(3,1,1)
plot(rf, Scheme(1).w_sc,'r-','DisplayName',sprintf('c=%0.1f,rf0=%0.4f,w2=%0.4f',Scheme(1).c_sc,Scheme(1).rf0_sc,Scheme(1).w2_sc));
hold on
for i=2:1:length(c)
color=[0.1*i,1-0.1*i,1-0.1*i];
plot(rf, Scheme(i).w_sc,'color', color,'DisplayName',sprintf('c=%0.1f,rf0=%0.4f,w2=%0.4f',Scheme(i).c_sc,Scheme(i).rf0_sc,Scheme(i).w2_sc));
end
xlim([0, 1]); 
ylim([0.6, 1.2]);
xlabel('rf');
ylabel('w');
title('a)rf-w曲线');
legend
subplot(3,1,2)
plot(rf, Scheme(1).Q_sc,'r-','DisplayName',sprintf('c=%0.1f,rf0=%0.4f,w2=%0.4f',Scheme(1).c_sc,Scheme(1).rf0_sc,Scheme(1).w2_sc));
hold on
for i=2:1:length(c)
color=[0.1*i,1-0.1*i,1-0.1*i];
plot(rf, Scheme(i).Q_sc,'color', color,'DisplayName',sprintf('c=%0.1f,rf0=%0.4f,w2=%0.4f',Scheme(i).c_sc,Scheme(i).rf0_sc,Scheme(i).w2_sc));
end
xlim([0, 1]); 
ylim([0, 200]);
xlabel('rf');
ylabel('Q');
title('b)rf-Q曲线');
legend
subplot(3,1,3)
plot(rf, Scheme(1).q_sc,'r-','DisplayName',sprintf('c=%0.1f,rf0=%0.4f,w2=%0.4f',Scheme(1).c_sc,Scheme(1).rf0_sc,Scheme(1).w2_sc));
hold on
for i=2:1:length(c)
color=[0.1*i,1-0.1*i,1-0.1*i];
plot(rf, Scheme(i).q_sc,'color', color,'DisplayName',sprintf('c=%0.1f,rf0=%0.4f,w2=%0.4f',Scheme(i).c_sc,Scheme(i).rf0_sc,Scheme(i).w2_sc));
end
xlim([0, 1]); 
ylim([0, 50]);
xlabel('rf');
ylabel('q');
title('c)rf-q曲线');
legend
%%均衡利润曲线
figure(4)
set(gcf, 'Position', [480, 100, 600, 600]);%%图像居中
subplot(4,1,1)
plot(rf, Scheme(1).S_sc,'r-','DisplayName',sprintf('c=%0.1f,rf0=%0.4f,w2=%0.4f',Scheme(1).c_sc,Scheme(1).rf0_sc,Scheme(1).w2_sc));
hold on
for i=2:1:length(c)
color=[0.1*i,1-0.1*i,1-0.1*i];
plot(rf, Scheme(i).S_sc,'color', color,'DisplayName',sprintf('c=%0.1f,rf0=%0.4f,w2=%0.4f',Scheme(i).c_sc,Scheme(i).rf0_sc,Scheme(i).w2_sc));
end
xlim([0, 1]); 
ylim([0, 60]);
xlabel('rf');
ylabel('S');
title('a)rf-S曲线');
legend
subplot(4,1,2)
plot(rf, Scheme(1).M_sc,'r-','DisplayName',sprintf('c=%0.1f,rf0=%0.4f,w2=%0.4f',Scheme(1).c_sc,Scheme(1).rf0_sc,Scheme(1).w2_sc));
hold on
for i=2:1:length(c)
color=[0.1*i,1-0.1*i,1-0.1*i];
plot(rf, Scheme(i).M_sc,'color', color,'DisplayName',sprintf('c=%0.1f,rf0=%0.4f,w2=%0.4f',Scheme(i).c_sc,Scheme(i).rf0_sc,Scheme(i).w2_sc));
end
xlim([0, 1]); 
ylim([0, 100]);
xlabel('rf');
ylabel('M');
title('b)rf-M曲线');
legend
subplot(4,1,3)
plot(rf, Scheme(1).R_sc,'r-','DisplayName',sprintf('c=%0.1f,rf0=%0.4f,w2=%0.4f',Scheme(1).c_sc,Scheme(1).rf0_sc,Scheme(1).w2_sc));
hold on
for i=2:1:length(c)
color=[0.1*i,1-0.1*i,1-0.1*i];
plot(rf, Scheme(i).R_sc,'color', color,'DisplayName',sprintf('c=%0.1f,rf0=%0.4f,w2=%0.4f',Scheme(i).c_sc,Scheme(i).rf0_sc,Scheme(i).w2_sc));
end
xlim([0, 1]); 
ylim([-1, 4]);
xlabel('rf');
ylabel('R');
title('c)rf-R曲线');
legend
subplot(4,1,4)
plot(rf, Scheme(1).L_sc,'r-','DisplayName',sprintf('c=%0.1f,rf0=%0.4f,w2=%0.4f',Scheme(1).c_sc,Scheme(1).rf0_sc,Scheme(1).w2_sc));
hold on
for i=2:1:length(c)
    color=[0.1*i,1-0.1*i,1-0.1*i];
plot(rf, Scheme(i).L_sc,'color', color,'DisplayName',sprintf('c=%0.1f,rf0=%0.4f,w2=%0.4f',Scheme(i).c_sc,Scheme(i).rf0_sc,Scheme(i).w2_sc));
end
xlim([0, 1]); 
ylim([40, 100]);
xlabel('rf');
ylabel('L');
title('c)rf-L曲线');
legend


%%
