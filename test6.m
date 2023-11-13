%%清理工作区数据
clear;                %清除所有变量
close all;                %清图
clc;
% 定义常数 a 和 k
a=100;
k=2;
x=1;
time = 0:0.1:1;
fx =(k/a)*(time/a).^(k-1).*exp(-(time/a).^k);%%韦伯分布密度
ft=1-exp(-(time./a).^k);%%韦伯分布
fk=exp(-(time./a).^k);%%韦伯反分布
fp=a*(-log(time)).^(1/k);%%韦伯反逆分布
subplot(411);plot(time,fx,'linewidth',2.5);xlabel('Time');ylabel('韦伯分布密度函数');
subplot(412);plot(time,ft,'linewidth',2.5);xlabel('Time');ylabel('韦伯分布函数');
subplot(413);plot(time,fk,'linewidth',2.5);xlabel('Time');ylabel('韦伯反分布');
subplot(414);plot(time,fp,'linewidth',2.5);xlabel('Time');ylabel('韦伯反逆分布');