function [rf0]=F_rf0_c(c,a,k)
equation = @(rf0) F_fan((F_fan(H_inverse(rf0/(1+rf0),a,k),a,k)/c*H_inverse(rf0/(1+rf0),a,k)),a,k)/c-1-rf0;
% �ṩ��ʼ�²�ֵ
initial_guess = 0.5; % �����ṩһ�����ʵĳ�ʼ�²�ֵ
% ʹ�� fsolve ��ⷽ��
rf0= fsolve(equation, initial_guess);