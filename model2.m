%%model2
%%������������
clear;                %������б���
close all;                %��ͼ
clc;
%%��������
%%��Ӧ��
%%������
%%������
%%Τ���ֲ���������
a=100;%%��������
k=2;%%��״���� 
%%���������������
c=0.2:0.05:0.8;
rf0_solution=zeros(1,length(c));
w1_solution=zeros(1,length(c));
for i=1:1:length(c)
    cp=c(i);
F=@(rf0,w1) cp*(1+rf0)*(F_fan_inverse(cp,a,k)-F_fan_inverse(cp*(1+rf0)/w1,a,k))+(1-w1)*integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+w1*integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k))-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp,a,k));
G=@(rf0,w1) (1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k));
% ��ʼ�²�ֵ
initial_guess = [1.0, 1.0];
equations = @(x) [F(x(1), x(2)); G(x(1), x(2))];
%options = optimoptions('fsolve', 'Display', 'iter');  % ������ʾ��������
solution = fsolve(equations,initial_guess);
rf0_solution(i)= solution(1);
w1_solution(i) = solution(2);
end
figure(1)
xlim([min(c), max(c)]); 
plot(c,rf0_solution,'r*-','DisplayName','c-rf0����');
hold on
xlim([min(c), max(c)]); 
plot(c,w1_solution,'bo-','DisplayName','c-w1����');
legend
clc;
%�������
%��������
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
%%��ͼ
figure(2)
set(gcf, 'Position', [480, 100, 600, 600]);%%ͼ�����
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
title('a)rf-w����');
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
title('b)rf-Q����');
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
title('c)rf-K����');
legend
%%������������
figure(3)
set(gcf, 'Position', [480, 100, 600, 600]);%%ͼ�����
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
title('a)rf-S����');
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
title('b)rf-M����');
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
title('c)rf-L����');
legend


