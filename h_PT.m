function [hp]=h_PT(x,a,k)
hp=(k/a)*(x/a).^(k-1);