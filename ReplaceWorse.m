%%最差个体统计并替代
function [chrom_new, fitness_new] = ReplaceWorse(chrom, chrom_best, fitness)

max_num = max(fitness);%%当代最高适应度
min_num = min(fitness);%%当代最低适应度

limit = (max_num-min_num)*0.8+min_num;%%设定适应度筛选阈值

replace_corr = fitness>limit;%%逻辑判断大于阈值为1小于阈值为0

replace_num = sum(replace_corr);%%统计大于阈值的个体个数
chrom(replace_corr,:) = ones(replace_num, 1)*chrom_best(1:2);%%将不合格个体染色体换成当代最优染色体
fitness(replace_corr) = ones(replace_num, 1)*chrom_best(end);%%将不合格个体适应度换成当代最优适应度
chrom_new = chrom;%%将筛选后的基因组赋给下一代的初始种群
fitness_new = fitness;
end