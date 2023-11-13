%%������ͳ�Ʋ����
function [chrom_new, fitness_new] = ReplaceWorse(chrom, chrom_best, fitness)

max_num = max(fitness);%%���������Ӧ��
min_num = min(fitness);%%���������Ӧ��

limit = (max_num-min_num)*0.8+min_num;%%�趨��Ӧ��ɸѡ��ֵ

replace_corr = fitness>limit;%%�߼��жϴ�����ֵΪ1С����ֵΪ0

replace_num = sum(replace_corr);%%ͳ�ƴ�����ֵ�ĸ������
chrom(replace_corr,:) = ones(replace_num, 1)*chrom_best(1:2);%%�����ϸ����Ⱦɫ�廻�ɵ�������Ⱦɫ��
fitness(replace_corr) = ones(replace_num, 1)*chrom_best(end);%%�����ϸ������Ӧ�Ȼ��ɵ���������Ӧ��
chrom_new = chrom;%%��ɸѡ��Ļ����鸳����һ���ĳ�ʼ��Ⱥ
fitness_new = fitness;
end