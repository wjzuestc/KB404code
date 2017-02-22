% @Description EDM算法
% @Author YD Wang
% @time 2016.12.22
%%
clc;
clear all;
close all;
%% 读入数据
% SampleData(:,1) = xlsread('coal','A:A');   %时间
SampleData(:,2) = xlsread('coal','B:B');   %数值
%%
figure;
plot(SampleData(:,2));