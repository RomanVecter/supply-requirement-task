function [cp,rf0,w1,RFMAX,WT,QT,KT,ST,MT,LT]=Model2_FUC(c,a,k)
cp=c;
rf0=0;
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


%%
rf=0:0.01:RFMAX;
WT=zeros(1,length(rf));
QT=zeros(1,length(rf));
KT=zeros(1,length(rf));
ST=zeros(1,length(rf));
MT=zeros(1,length(rf));
LT=zeros(1,length(rf));
for i=1:1:(length(rf))
    rdp=rf(i);
    %%���ݵ�ǰ��c��rf0�����w1
    equation = @(w1) (1/w1-1)*c/w1/h_PT(F_fan_inverse(c/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(c/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(c*(1+rdp)/w1,a,k));
    % �ṩ��ʼ�²�ֵ
    initial_guess = 0.81; % �����ṩһ�����ʵĳ�ʼ�²�ֵ
    % ʹ�� fsolve ��ⷽ��
    w1= fsolve(equation, initial_guess);
    w=w1;
    t=c/w;
    p=c*(1+rdp);
    e=c*(1+rdp)/w;
    if p>=1
        fprintf('���ڴ�������ʹ��c*(1+rf)>=1,��ʱc=%0.1f,rf0=%0.4f,w1=%0.4f,rf=%0.3f\n',c,rf0,w1,rdp);
        %disp(['��' sprintf( '%i',i) '��������']);
    end
    if rdp<=rf0
       %%�����
       WT(i)=w;
       QT(i)=F_fan_inverse(c,a,k);
       KT(i)=F_fan_inverse(c,a,k);
       %%����������
       ST(i)=0;
       %%M���ֲ�������
       fun1=@(x) F_fan(x,a,k)-p;  %%���ֺ���1
       a1=0; %%��������
       b1=F_fan_inverse(c,a,k);%%��������
       MT(i)=integral(fun1, a1, b1);
       %%L���ֲ�������
       fun2=@(x) F_fan(x,a,k)-p;  %%���ֺ���1
       a2=0; %%��������
       b2=F_fan_inverse(c,a,k);%%��������
       LT(i)=integral(fun2,a2, b2);
    else
       %%�����
       WT(i)=w;
       QT(i)=F_fan_inverse(e,a,k);
       KT(i)=F_fan_inverse(t,a,k);
       %%����������
       %%S���ֲ�������
       fun1=@(x) F_distribution(x,a,k);  %%���ֺ���1
       a1=F_fan_inverse(e,a,k); %%��������
       b1=F_fan_inverse(t,a,k);%%��������
       ST(i)=(w-c)*(F_fan_inverse(t,a,k)-F_fan_inverse(e,a,k))-w*integral(fun1,a1, b1);
       %%M���ֲ�������
       fun1=@(x) F_fan(x,a,k);  %%���ֺ���1
       a1=0; %%��������
       b1=F_fan_inverse(t,a,k);%%��������
       fun2=@(x) w*F_fan(x,a,k)-p;%%���ֺ���2
       a2=0; %%��������
       b2=F_fan_inverse(e,a,k);%%��������
       MT(i)=(1-w)*integral(fun1,a1, b1)+integral(fun2,a2, b2);
       %%L���ֲ�������
       fun1=@(x) F_fan(x,a,k);  %%���ֺ���1
       a1=0; %%��������
       b1=F_fan_inverse(t,a,k);%%��������
       LT(i)=integral(fun1,a1, b1)-c*(F_fan_inverse(t,a,k)+rdp*F_fan_inverse(e,a,k));
    end
end
%%��w1_min
%%���ݵ�ǰ��c��RFMAX�����w1_min
        equation = @(w1) (1/w1-1)*c/w1/h_PT(F_fan_inverse(c/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(c/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(c*(1+RFMAX)/w1,a,k));
        % �ṩ��ʼ�²�ֵ
        initial_guess = 0.81; % �����ṩһ�����ʵĳ�ʼ�²�ֵ
        % ʹ�� fsolve ��ⷽ��
        w1= fsolve(equation, initial_guess);
