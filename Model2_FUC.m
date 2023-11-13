function [cp,rf0,w1,RFMAX,WT,QT,KT,ST,MT,LT]=Model2_FUC(c,a,k)
cp=c;
rf0=0;
%%
RF=0:0.01:1;
W=zeros(1,length(RF));
CP=zeros(1,length(RF));
for j=1:1:length(RF)
    rdp=RF(j);
% 定义函数
equation = @(w1) (1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rdp)/w1,a,k));
% 提供初始猜测值
initial_guess = 0.81; % 可以提供一个合适的初始猜测值

% 使用 fsolve 求解方程
W(j)= fsolve(equation, initial_guess);
if (j~=1)&&(W(j)>W(j-1))
 initial_guess = 0.1; % 可以提供一个合适的初始猜测值
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
    %%根据当前的c和rf0，求解w1
    equation = @(w1) (1/w1-1)*c/w1/h_PT(F_fan_inverse(c/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(c/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(c*(1+rdp)/w1,a,k));
    % 提供初始猜测值
    initial_guess = 0.81; % 可以提供一个合适的初始猜测值
    % 使用 fsolve 求解方程
    w1= fsolve(equation, initial_guess);
    w=w1;
    t=c/w;
    p=c*(1+rdp);
    e=c*(1+rdp)/w;
    if p>=1
        fprintf('存在错误输入使得c*(1+rf)>=1,此时c=%0.1f,rf0=%0.4f,w1=%0.4f,rf=%0.3f\n',c,rf0,w1,rdp);
        %disp(['第' sprintf( '%i',i) '发生干涉']);
    end
    if rdp<=rf0
       %%均衡解
       WT(i)=w;
       QT(i)=F_fan_inverse(c,a,k);
       KT(i)=F_fan_inverse(c,a,k);
       %%均衡利润函数
       ST(i)=0;
       %%M积分参数设置
       fun1=@(x) F_fan(x,a,k)-p;  %%积分函数1
       a1=0; %%积分下限
       b1=F_fan_inverse(c,a,k);%%积分上限
       MT(i)=integral(fun1, a1, b1);
       %%L积分参数设置
       fun2=@(x) F_fan(x,a,k)-p;  %%积分函数1
       a2=0; %%积分下限
       b2=F_fan_inverse(c,a,k);%%积分上限
       LT(i)=integral(fun2,a2, b2);
    else
       %%均衡解
       WT(i)=w;
       QT(i)=F_fan_inverse(e,a,k);
       KT(i)=F_fan_inverse(t,a,k);
       %%均衡利润函数
       %%S积分参数设置
       fun1=@(x) F_distribution(x,a,k);  %%积分函数1
       a1=F_fan_inverse(e,a,k); %%积分下限
       b1=F_fan_inverse(t,a,k);%%积分上限
       ST(i)=(w-c)*(F_fan_inverse(t,a,k)-F_fan_inverse(e,a,k))-w*integral(fun1,a1, b1);
       %%M积分参数设置
       fun1=@(x) F_fan(x,a,k);  %%积分函数1
       a1=0; %%积分下限
       b1=F_fan_inverse(t,a,k);%%积分上限
       fun2=@(x) w*F_fan(x,a,k)-p;%%积分函数2
       a2=0; %%积分下限
       b2=F_fan_inverse(e,a,k);%%积分上限
       MT(i)=(1-w)*integral(fun1,a1, b1)+integral(fun2,a2, b2);
       %%L积分参数设置
       fun1=@(x) F_fan(x,a,k);  %%积分函数1
       a1=0; %%积分下限
       b1=F_fan_inverse(t,a,k);%%积分上限
       LT(i)=integral(fun1,a1, b1)-c*(F_fan_inverse(t,a,k)+rdp*F_fan_inverse(e,a,k));
    end
end
%%求w1_min
%%根据当前的c和RFMAX，求解w1_min
        equation = @(w1) (1/w1-1)*c/w1/h_PT(F_fan_inverse(c/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(c/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(c*(1+RFMAX)/w1,a,k));
        % 提供初始猜测值
        initial_guess = 0.81; % 可以提供一个合适的初始猜测值
        % 使用 fsolve 求解方程
        w1= fsolve(equation, initial_guess);
