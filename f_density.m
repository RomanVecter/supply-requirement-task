%%韦伯分布概率密度函数
function [f]=f_density(x,a,k)
f=(k/a)*(x/a)^(k-1)*exp(-(x/a)^k);