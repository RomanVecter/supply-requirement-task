%%model5
%%��֤��Ӧ����������������������ͼ
%%�޸�15�е�c������ʵ�ֲ�ͬc��Ӧ�������ͼ��
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
%%�Լ��趨c
c=0.3;

[cp,rf0,w2,w1_min,RFMAX,L_F,M_F,R_F,L_B,M_B,R_B]=Model5_FUC(c,a,k);
L_limit=0;
H_limit=100;

%%
x1=rf0;
x2=RFMAX;
x_limit=x2;
figure(1)
set(gcf, 'Position', [480, 100, 1000, 600]);%%ͼ�����
rf=0.01:0.01:1;
plot(x1* ones(size(rf)), rf*100,'r--','DisplayName',sprintf('rf=rf0,    c=%0.1f,    rf0=%0.2f',c,rf0),'LineWidth', 2 );  % ���ƺ�ɫֱ��
hold on
plot(x2* ones(size(rf)), rf*100, 'b--','DisplayName',sprintf('rf=RFMAX,     c=%0.1f,    RFMAX=%0.2f',c,x2),'LineWidth', 2);  % ������ɫֱ��
xlim([0, x_limit*1.2]); 
rf=0.01:0.01:RFMAX;
plot(rf, L_F,'color', [0.1,0.9,0.5],'DisplayName',sprintf('��Ӧ����������L_F,c=%0.1f*(1+RFMAX=%0.2f)<w1=%0.2f',c,RFMAX,w1_min),'LineWidth',2); 
plot(rf, M_F,'color', [0.5,0.5,0.3],'DisplayName',sprintf('��Ӧ����������M_F,c=%0.1f*(1+RFMAX=%0.2f)<w1=%0.2f',c,RFMAX,w1_min),'LineWidth',2);
plot(rf, R_F,'color', [0.8,0.6,0.2],'DisplayName',sprintf('��Ӧ����������R_F,c=%0.1f*(1+RFMAX=%0.2f)<w1=%0.2f',c,RFMAX,w1_min),'LineWidth',2); 

xlabel('rf');
ylabel('��Ӧ������');
title('��Ӧ����������');
legend

figure(2)
set(gcf, 'Position', [480, 100, 1000, 600]);%%ͼ�����
rf=0.01:0.01:1;
plot(x1* ones(size(rf)), rf*100,'r--','DisplayName',sprintf('rf=rf0,    c=%0.1f,    rf0=%0.2f',c,rf0),'LineWidth', 2 );  % ���ƺ�ɫֱ��
hold on
plot(x2* ones(size(rf)), rf*100, 'b--','DisplayName',sprintf('rf=RFMAX,     c=%0.1f,    RFMAX=%0.2f',c,x2),'LineWidth', 2);  % ������ɫֱ��
xlim([0, x_limit*1.2]); 
rf=0.01:0.01:RFMAX;
plot(rf, L_B,'color', [0.1,0.9,0.5],'DisplayName',sprintf('��������������L_B,c=%0.1f*(1+RFMAX=%0.2f)<w1=%0.2f',c,RFMAX,w1_min),'LineWidth',2);  
plot(rf, M_B,'color', [0.5,0.5,0.3],'DisplayName',sprintf('��������������M_B,c=%0.1f*(1+RFMAX=%0.2f)<w1=%0.2f',c,RFMAX,w1_min),'LineWidth',2);
plot(rf, R_B,'color', [0.8,0.6,0.2],'DisplayName',sprintf('��������������R_B,c=%0.1f*(1+RFMAX=%0.2f)<w1=%0.2f',c,RFMAX,w1_min),'LineWidth',2); 

xlabel('rf');
ylabel('����������');
title('��������������');
legend