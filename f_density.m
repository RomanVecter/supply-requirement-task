%%Τ���ֲ������ܶȺ���
function [f]=f_density(x,a,k)
f=(k/a)*(x/a)^(k-1)*exp(-(x/a)^k);