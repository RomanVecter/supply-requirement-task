%%������������
clear;                %������б���
close all;                %��ͼ
clc;
% ���峣�� a �� k
a=100;
k=2;
cp=0.2;
rf0=0;
w1=1;
P=(1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k));
cp*(1+rf0)/w1

F_fan_inverse(cp*(1+rf0)/w1,a,k)
(1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)
integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))
integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k))