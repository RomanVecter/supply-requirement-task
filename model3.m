%%model3
%%���һЩ��Ҫ��֤����
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
%%��֤1
CP=0.2:0.01:0.6;
WP=zeros(1,length(CP));
for i=1:1:length(CP)
    cp=CP(i);
equation = @(w0) (1/w0-1)*cp/w0/h_PT(F_fan_inverse(cp/w0,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w0,a,k));
% �ṩ��ʼ�²�ֵ
initial_guess = 0.95; % �����ṩһ�����ʵĳ�ʼ�²�ֵ
% ʹ�� fsolve ��ⷽ��
WP(i)= fsolve(equation, initial_guess);
end
figure(1)
plot(CP,WP,'r*-','DisplayName','c-w0����');
hold on
plot(CP,CP,'bo-','DisplayName','w0=c');
xlabel('c');
ylabel('w0');
legend
%%��֤2
a = 100;
k = 2;
cp = 0.5;
equation = @(w0) (1/w0-1)*cp/w0/h_PT(F_fan_inverse(cp/w0,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w0,a,k));
% �ṩ��ʼ�²�ֵ
initial_guess = 0.1; % �����ṩһ�����ʵĳ�ʼ�²�ֵ
% ʹ�� fsolve ��ⷽ��
w0= fsolve(equation, initial_guess);
W=w0+0.01:0.01:1;
RF=zeros(1,length(W));
for i=1:1:length(W)
    w1=W(i);
% ���庯��
equation = @(rf0) (1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k));
% �ṩ��ʼ�²�ֵ
initial_guess = 0.1; % �����ṩһ�����ʵĳ�ʼ�²�ֵ
% ʹ�� fsolve ��ⷽ��
RF(i)= fsolve(equation, initial_guess);
end
HDP=(1/cp-1)*ones(1,length(RF));
figure(2)
plot(RF,W,'r*-','DisplayName','w1-rf����');
hold on
plot(HDP,W,'bo-','DisplayName','rf=1/cp-1');
xlim([0, 1]); 
% ylim([0, 1]);
xlabel('rf');
ylabel('w1');
legend

%%��֤3
a = 100;
k = 2;
RFP=WP./CP-1;
HDP=1./CP-1;
figure(3)
plot(CP,RFP,'r*-','DisplayName','c-rf����');
hold on
plot(CP,HDP,'bo-','DisplayName','1/cp-1');
xlabel('c');
ylabel('rf');
legend

WT=zeros(1,length(CP));
for i=1:1:length(CP)
    cp=CP(i);
    rf0=RFP(i);
equation = @(w1) (1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k));
% �ṩ��ʼ�²�ֵ
initial_guess = 0.8; % �����ṩһ�����ʵĳ�ʼ�²�ֵ
% ʹ�� fsolve ��ⷽ��
WT(i)= fsolve(equation, initial_guess);
end
figure(4)
plot(CP,WT,'r*-','DisplayName','c-w1����');
hold on
plot(CP,CP,'bo-','DisplayName','w1=c');
xlabel('c');
ylabel('w1');
legend

%%��֤4
a = 100;
k = 2;
ZO=zeros(1,length(CP));
for i=1:1:length(CP)
    cp=CP(i);
    rf0=RFP(i);
    w1=WT(i);
ZO(i)=cp*(1+rf0)*(F_fan_inverse(cp,a,k)-F_fan_inverse(cp*(1+rf0)/w1,a,k))+(1-w1)*integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+w1*integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k))-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp,a,k));
if ZO(i)<0
    fprintf('���ڴ�������ʹ��ʽ��С�ڵ���0,��ʱc=%0.4f,w1=%0.4f,rf=%0.4f\n',cp,w1,rf0);
end
end
figure(5)
ZT=zeros(1,length(CP));
plot(CP,ZO,'r*-','DisplayName','c-F����');
hold on
plot(CP,ZT,'bo-','DisplayName','F=0');
xlabel('c');
ylabel('F');
legend