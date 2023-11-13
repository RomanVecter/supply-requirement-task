%%������������
clear;                %������б���
close all;                %��ͼ
clc;
% ���峣�� a �� k
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
% ���庯��
equation = @(w1) (1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rdp)/w1,a,k));
% �ṩ��ʼ�²�ֵ
initial_guess =0.82; % �����ṩһ�����ʵĳ�ʼ�²�ֵ

% ʹ�� fsolve ��ⷽ��
W(j)= fsolve(equation, initial_guess);
if (j~=1)&&(W(j)>W(j-1))
 initial_guess = 0.1; % �����ṩһ�����ʵĳ�ʼ�²�ֵ
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