% @Description EDM�㷨
% @Author YD Wang
% @time 2016.12.22
%%
clc;
clear all;
close all;
%% ��������
% SampleData(:,1) = xlsread('coal','A:A');   %ʱ��
SampleData(:,2) = xlsread('coal','B:B');   %��ֵ
%%
figure;
plot(SampleData(:,2));