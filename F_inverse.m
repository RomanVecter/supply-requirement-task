%%Τ���ֲ�������
function [F_ni]=F_inverse(y,a,k)
F_ni=a*(-log(1-y))^(1/k);