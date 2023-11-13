function [cp,rf0,w2,w1_min,RFMAX,L_F,M_F,R_F,L_B,M_B,R_B]=Model5_FUC(c,a,k)
cp=c;
rf0=F_rf0_c(c,a,k);
w2=F_w2_c(c,a,k);
%%
RF=0:0.01:1;
W=zeros(1,length(RF));
CP=zeros(1,length(RF));
for j=1:1:length(RF)
    rdp=RF(j);
% ���庯��
equation = @(w1) (1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rdp)/w1,a,k));
% �ṩ��ʼ�²�ֵ
initial_guess = 0.81; % �����ṩһ�����ʵĳ�ʼ�²�ֵ

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

%%
%%��w1_min
%%���ݵ�ǰ��c��RFMAX�����w1_min
        equation = @(w1) (1/w1-1)*c/w1/h_PT(F_fan_inverse(c/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(c/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(c*(1+RFMAX)/w1,a,k));
        % �ṩ��ʼ�²�ֵ
        initial_guess = 0.81; % �����ṩһ�����ʵĳ�ʼ�²�ֵ
        % ʹ�� fsolve ��ⷽ��
        w1_min= fsolve(equation, initial_guess);
%%
rf=0.01:0.01:RFMAX;
L_F=zeros(1,length(rf));
M_F=zeros(1,length(rf));
R_F=zeros(1,length(rf));
L_B=zeros(1,length(rf));
M_B=zeros(1,length(rf));
R_B=zeros(1,length(rf));
    for i=1:1:(length(rf))
        rdp=rf(i);
        t=rdp/(1+rdp);
        p=c*(1+rdp);
        %%���ݵ�ǰ��c��rf0�����w1
        equation = @(w1) (1/w1-1)*c/w1/h_PT(F_fan_inverse(c/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(c/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(c*(1+rdp)/w1,a,k));
        % �ṩ��ʼ�²�ֵ
        initial_guess = 0.81; % �����ṩһ�����ʵĳ�ʼ�²�ֵ
        % ʹ�� fsolve ��ⷽ��
        w1= fsolve(equation, initial_guess);
        %%w1������
        r=c/w1;
        e=c*(1+rdp)/w1;
        if rdp<rf0
           %%��Ӧ������F 
        %%L���ֲ�������
        fun3=@(x) F_fan(x,a,k);%%���ֺ���2
        a3=0; %%��������
        b3=F_fan_inverse(p,a,k);%%��������
        L_F(i)=rdp*(F_fan(H_inverse(t,a,k),a,k)*H_inverse(t,a,k)-c*F_fan_inverse(p,a,k))+integral(fun3, a3, b3);
        %%M���ֲ�������
        fun1=@(x) F_fan(x,a,k)-c;  %%���ֺ���1
        a1=0; %%��������
        b1=F_fan_inverse(c,a,k);%%��������
        M_F(i)=integral(fun1, a1, b1);
        %%R���ֲ�������
        fun1=@(x) F_fan(x,a,k);  %%���ֺ���1
        a1=0; %%��������
        b1=F_fan_inverse(r,a,k);%%��������
        R_F(i)=integral(fun1,a1, b1)-c*(F_fan_inverse(r,a,k)+rdp*F_fan_inverse(e,a,k));
        
        %%������_B
        %%L���ֲ�������
        fun1=@(x) F_fan(x,a,k)-p;  %%���ֺ���1
        a1=0; %%��������
        b1=F_fan_inverse(p,a,k);%%��������
        fun2=@(x) F_fan(x,a,k)-F_fan(H_inverse(t,a,k),a,k);%%���ֺ���2
        a2=0; %%��������
        b2=H_inverse(t,a,k);%%��������
        L_B(i)=integral(fun1, a1, b1)-integral(fun2, a2, b2);
        %%M���ֲ�������
        fun1=@(x) F_fan(x,a,k)-c;  %%���ֺ���1
        a1=0; %%��������
        b1=F_fan_inverse(c,a,k);%%��������
        M_B(i)=integral(fun1, a1, b1);
        %%R���ֲ�������
       fun1=@(x) F_fan(x,a,k);  %%���ֺ���1
       a1=0; %%��������
       b1=F_fan_inverse(r,a,k);%%��������
       fun2=@(x) w1*F_fan(x,a,k)-p;%%���ֺ���2
       a2=0; %%��������
       b2=F_fan_inverse(e,a,k);%%��������
       R_B(i)=(1-w1)*integral(fun1,a1, b1)+integral(fun2,a2, b2);
        else
            
            %%��Ӧ������F
        %%L���ֲ�������
        fun3=@(x) F_fan(x,a,k);%%���ֺ���2
        a3=0; %%��������
        b3=w2/c*F_fan_inverse(w2,a,k);%%��������
        L_F(i)=integral(fun3, a3, b3);    
        %%M���ֲ�������
        fun1=@(x) F_fan(x,a,k)-c;  %%���ֺ���1
        a1=0; %%��������
        b1=F_fan_inverse(c,a,k);%%��������
        M_F(i)=integral(fun1, a1, b1);
        %%R���ֲ�������
        fun1=@(x) F_fan(x,a,k);  %%���ֺ���1
        a1=0; %%��������
        b1=F_fan_inverse(r,a,k);%%��������
        R_F(i)=integral(fun1,a1, b1)-c*(F_fan_inverse(r,a,k)+rdp*F_fan_inverse(e,a,k));
        
        %%������_B
        %%L���ֲ�������
        fun1=@(x) F_distribution(x,a,k);  %%���ֺ���1
        a1=F_fan_inverse(w2,a,k); %%��������
        b1=w2/c*F_fan_inverse(w2,a,k);%%��������
        L_B(i)=(w2/c-1)*F_fan_inverse(w2,a,k)-integral(fun1, a1, b1);
        %%M���ֲ�������
        fun1=@(x) F_fan(x,a,k)-c;  %%���ֺ���1
        a1=0; %%��������
        b1=F_fan_inverse(c,a,k);%%��������
        M_B(i)=integral(fun1, a1, b1);
        %%R���ֲ�������
       fun1=@(x) F_fan(x,a,k);  %%���ֺ���1
       a1=0; %%��������
       b1=F_fan_inverse(r,a,k);%%��������
       fun2=@(x) w1*F_fan(x,a,k)-p;%%���ֺ���2
       a2=0; %%��������
       b2=F_fan_inverse(e,a,k);%%��������
       R_B(i)=(1-w1)*integral(fun1,a1, b1)+integral(fun2,a2, b2);
        end
        
    end
