function [w2]=F_w2_c(c,a,k)
equation = @(w2) F_fan(w2/c*F_fan_inverse(w2,a,k),a,k)*(1-H(F_fan_inverse(w2,a,k),a,k))/c-1;
% �ṩ��ʼ�²�ֵ
initial_guess = 0.5; % �����ṩһ�����ʵĳ�ʼ�²�ֵ
% ʹ�� fsolve ��ⷽ��
w2= fsolve(equation, initial_guess);