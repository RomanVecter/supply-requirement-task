%%������������
clear;                %������б���
close all;                %��ͼ
clc;
% ���峣�� a �� k
a=100;
k=2;
x=1;
time = 0:0.1:1;
fx =(k/a)*(time/a).^(k-1).*exp(-(time/a).^k);%%Τ���ֲ��ܶ�
ft=1-exp(-(time./a).^k);%%Τ���ֲ�
fk=exp(-(time./a).^k);%%Τ�����ֲ�
fp=a*(-log(time)).^(1/k);%%Τ������ֲ�
subplot(411);plot(time,fx,'linewidth',2.5);xlabel('Time');ylabel('Τ���ֲ��ܶȺ���');
subplot(412);plot(time,ft,'linewidth',2.5);xlabel('Time');ylabel('Τ���ֲ�����');
subplot(413);plot(time,fk,'linewidth',2.5);xlabel('Time');ylabel('Τ�����ֲ�');
subplot(414);plot(time,fp,'linewidth',2.5);xlabel('Time');ylabel('Τ������ֲ�');