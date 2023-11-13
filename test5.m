%%������������
clear;                %������б���
close all;                %��ͼ
clc;
% ���峣�� a �� k
a=100;
k=2;
cp=0.5;
rf0=0;
equation = @(w1) (1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k));
% �ṩ��ʼ�²�ֵ
initial_guess =1; % �����ṩһ�����ʵĳ�ʼ�²�ֵ
% ʹ�� fsolve ��ⷽ��

w1= fsolve(equation, initial_guess);
if isreal(w1)
    disp('w1����һ��ʵ����');
else
    disp('w1����һ��������');
end
cp*(1+rf0)/w1

F_fan_inverse(cp*(1+rf0)/w1,a,k)
(1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)
integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))
integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k))
FP=(1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k))