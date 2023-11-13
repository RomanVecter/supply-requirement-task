%%寻找当代最优适应度
function chrom_best= FindBest(chrom,fitness, N_chrom,cp,a,k)
chrom_best = zeros(1, N_chrom+3);
[maxNum, maxCorr] = min(fitness);%寻找最小适应度的大小及下标值
chrom_best(1:N_chrom) =chrom(maxCorr,:);%%存储当代最佳染色体
w1=chrom_best(1);
rf0=chrom_best(2);
chrom_best(3)=cp*(1+rf0)*(F_fan_inverse(cp,a,k)-F_fan_inverse(cp*(1+rf0)/w1,a,k))+(1-w1)*integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+w1*integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k))-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp,a,k));
chrom_best(4)=(1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k));    
chrom_best(end) = maxNum;%%存储当代最优适应度