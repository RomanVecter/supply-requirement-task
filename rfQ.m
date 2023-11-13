%%清理工作区数据
clear;                %清除所有变量
close all;                %清图
clc;
a=100;
k=2;
c=0.5;
global rf
i=1;
for rf=0:0.01:2
Q1(i)=F_inverse(1-c*(1+rf),a,k);
H_ni(i)=H_inverse(rf/(1+rf),a,k);
Q2(i)=(1-F_distribution(H_ni(i),a,k))*H_ni(i)/c;
i=i+1;
end
rf=0:0.01:2;
plot(rf,Q1,'color','r')
hold on
plot(rf,Q2)
