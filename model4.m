%%model4
%%利用遗传算法求解隐函数方程组，得到给定c下，解得w1与rf

%%清理工作区数据
clear;                %清除所有变量
close all;                %清图
clc;
%%参数定义
%%供应商
%%制造商
%%零售商
%%韦伯分布参数控制
a=100;%%缩放因子
k=2;%%形状参数 
cp=0.8;
%%遗传算法求解
%% 基础参数
N = 20;  %种群内个体数目
N_chrom = 2; %染色体节点数，也就是每个个体有多少条染色体，适应函数里有几个自变量。
iter = 150; %迭代次数，也就是一共有多少代
mut = 0.2;  %突变概率
acr = 0.2; %交叉概率
chrom_range = [cp 1;0 1];%每个节点的值的区间
%%参数存储
chrom = zeros(N, N_chrom);%族群染色体存储矩阵
F=zeros(N, 1);
G=zeros(N, 1);
fitness = zeros(N, 1);%族群每个个体当代染色体的适应度存储矩阵
fitness_ave = zeros(1, iter);%存放每一代的平均适应度
fitness_best = zeros(1, iter);%存放每一代的最优适应度
chrom_best = zeros(1, N_chrom+3);%存放当前代的最优染色体与适应度
%% 初始化，这只是用于生成第一代个体，并计算其适应度函数
for i=1:N
    for j=1:N_chrom 
        chrom(i,j)=unifrnd(chrom_range(j,1),chrom_range(j,2));
    end
end
[F,G,fitness]= CalFitness(chrom,N,cp,a,k);%%%计算每一个个体适应度
%%寻找初代最优染色体
chrom_best = FindBest(chrom,fitness, N_chrom,cp,a,k);
fitness_best(1) = chrom_best(end); %将当前最优适应度存入矩阵当中
[N ,~] = size(fitness);%%获取数组fitness的行数并赋给N
fitness_ave(1) = sum(fitness)/N;
for t = 2:iter
     chrom = MutChrom(chrom, mut, N, N_chrom, chrom_range, t, iter); %变异
     chrom = AcrChrom(chrom, acr, N, N_chrom); %交叉
     [F,G,fitness]= CalFitness(chrom,N,cp,a,k);%%%计算每一个个体适应度
     chrom_best_temp= FindBest(chrom,fitness, N_chrom,cp,a,k);%寻找最优染色体
    if chrom_best_temp(end)<chrom_best(end) %如果当代的最优适应度小于上一代储存的最优适应度
        chrom_best = chrom_best_temp;
    end
    [chrom, fitness] = ReplaceWorse(chrom, chrom_best, fitness);
    fitness_best(t) = chrom_best(end); %将当前最优适应度存入矩阵当中
    [N ,~] = size(fitness);%%获取数组fitness的行数并赋给N
    fitness_ave(t) = sum(fitness)/N;
end
%% 作图
figure(1)
plot(1:iter, fitness_best, 'b')
grid on
%% 输出结果
disp(['w1和rf分别为：', num2str(chrom_best(1:2))])
disp(['方程组计算误差最小为：', num2str(chrom_best(end))])