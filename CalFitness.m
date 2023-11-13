function [F,G,fitness]= CalFitness(chrom,N,cp,a,k)
F=zeros(N, 1);
G=zeros(N, 1);
fitness = zeros(N, 1);%族群每个个体当代染色体的适应度存储矩阵
for i=1:N
    w1=chrom(i,1);
    rf0=chrom(i,2);
    F(i)=cp*(1+rf0)*(F_fan_inverse(cp,a,k)-F_fan_inverse(cp*(1+rf0)/w1,a,k))+(1-w1)*integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+w1*integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k))-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp,a,k));
    G(i)=(1/w1-1)*cp/w1/h_PT(F_fan_inverse(cp/w1,a,k),a,k)-integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp/w1,a,k))+integral(@(x) F_fan(x,a,k),0,F_fan_inverse(cp*(1+rf0)/w1,a,k));    
    fitness(i)=abs(F(i))+abs(G(i));
    if fitness(i)>1000
        fitness(i)=1000;
    end
end

    
end