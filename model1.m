%%model1
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
%%������
syms c rf0
eqn=F_fan((F_fan(H_inverse(rf0/(1+rf0),a,k),a,k)/c*H_inverse(rf0/(1+rf0),a,k)),a,k)/c-1-rf0;
fimplicit(eqn, [0, 1, 0, 1],'DisplayName','rf0'); % ָ�� x �� y �ķ�Χ
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
%%c-w2����
syms c w2
epn=F_fan(w2/c*F_fan_inverse(w2,a,k),a,k)*(1-H(F_fan_inverse(w2,a,k),a,k))/c-1;
figure(2);
fimplicit(epn, [0, 1, 0, 1],'DisplayName','w2'); % ָ�� x �� y �ķ�Χ
hold on
c=0:0.01:1;
w2=c;
plot(c,w2,'DisplayName','w2=c');
xlabel('c');
ylabel('w2');
legend;

%%��������
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
%%��ͼ
rf=0.01:0.01:1;
figure(3)
set(gcf, 'Position', [480, 100, 600, 600]);%%ͼ�����
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
title('a)rf-w����');
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
title('b)rf-Q����');
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
title('c)rf-q����');
legend
%%������������
figure(4)
set(gcf, 'Position', [480, 100, 600, 600]);%%ͼ�����
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
title('a)rf-S����');
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
title('b)rf-M����');
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
title('c)rf-R����');
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
title('c)rf-L����');
legend


%%
