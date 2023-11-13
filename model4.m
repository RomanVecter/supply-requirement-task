%%model4
%%�����Ŵ��㷨��������������飬�õ�����c�£����w1��rf

%%������������
clear;                %������б���
close all;                %��ͼ
clc;
%%��������
%%��Ӧ��
%%������
%%������
%%Τ���ֲ���������
a=100;%%��������
k=2;%%��״���� 
cp=0.8;
%%�Ŵ��㷨���
%% ��������
N = 20;  %��Ⱥ�ڸ�����Ŀ
N_chrom = 2; %Ⱦɫ��ڵ�����Ҳ����ÿ�������ж�����Ⱦɫ�壬��Ӧ�������м����Ա�����
iter = 150; %����������Ҳ����һ���ж��ٴ�
mut = 0.2;  %ͻ�����
acr = 0.2; %�������
chrom_range = [cp 1;0 1];%ÿ���ڵ��ֵ������
%%�����洢
chrom = zeros(N, N_chrom);%��ȺȾɫ��洢����
F=zeros(N, 1);
G=zeros(N, 1);
fitness = zeros(N, 1);%��Ⱥÿ�����嵱��Ⱦɫ�����Ӧ�ȴ洢����
fitness_ave = zeros(1, iter);%���ÿһ����ƽ����Ӧ��
fitness_best = zeros(1, iter);%���ÿһ����������Ӧ��
chrom_best = zeros(1, N_chrom+3);%��ŵ�ǰ��������Ⱦɫ������Ӧ��
%% ��ʼ������ֻ���������ɵ�һ�����壬����������Ӧ�Ⱥ���
for i=1:N
    for j=1:N_chrom 
        chrom(i,j)=unifrnd(chrom_range(j,1),chrom_range(j,2));
    end
end
[F,G,fitness]= CalFitness(chrom,N,cp,a,k);%%%����ÿһ��������Ӧ��
%%Ѱ�ҳ�������Ⱦɫ��
chrom_best = FindBest(chrom,fitness, N_chrom,cp,a,k);
fitness_best(1) = chrom_best(end); %����ǰ������Ӧ�ȴ��������
[N ,~] = size(fitness);%%��ȡ����fitness������������N
fitness_ave(1) = sum(fitness)/N;
for t = 2:iter
     chrom = MutChrom(chrom, mut, N, N_chrom, chrom_range, t, iter); %����
     chrom = AcrChrom(chrom, acr, N, N_chrom); %����
     [F,G,fitness]= CalFitness(chrom,N,cp,a,k);%%%����ÿһ��������Ӧ��
     chrom_best_temp= FindBest(chrom,fitness, N_chrom,cp,a,k);%Ѱ������Ⱦɫ��
    if chrom_best_temp(end)<chrom_best(end) %���������������Ӧ��С����һ�������������Ӧ��
        chrom_best = chrom_best_temp;
    end
    [chrom, fitness] = ReplaceWorse(chrom, chrom_best, fitness);
    fitness_best(t) = chrom_best(end); %����ǰ������Ӧ�ȴ��������
    [N ,~] = size(fitness);%%��ȡ����fitness������������N
    fitness_ave(t) = sum(fitness)/N;
end
%% ��ͼ
figure(1)
plot(1:iter, fitness_best, 'b')
grid on
%% ������
disp(['w1��rf�ֱ�Ϊ��', num2str(chrom_best(1:2))])
disp(['��������������СΪ��', num2str(chrom_best(end))])