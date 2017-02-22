% @用前向欧拉的方式求解ode
% @Author JingzengWang
% @time 2016.11.15
%%
clc;
clear all;
close all;
%% 定义构造函数参量
func1 = @(x)sin(2*pi*x)+cos(2*pi*x);   
func2 = @(x)3*x+1;
func3 = @(x)x.*x+1;
f1_init = func1(0);           %设置初始点
f2_init = func2(0);
f3_init = func3(0);
constant_A = 2;               %微分方程常系数
sampling_step = 0.0001;         %采样间隔
% sampling_step = 1;          %采样间隔
% samp_max = 100;
samp_max = 10;
samp_min = 0;
samp_num = (samp_max-samp_min)/sampling_step;
sampling_x = linspace(samp_min,samp_max,samp_num);
%% 通过构造函数获得恒解空间
% func1_diff = diff(func1,'x');
% func2_diff = diff(func2,'x');
% func3_diff = diff(func3,'x');
func1_diff = @(x)2*pi*cos(2*pi*x)-2*pi*sin(2*pi*x);
func2_diff = @(x)3;
func3_diff = @(x)2*x;
for i=1:samp_num
	g1_fun(i) = func1(sampling_x(i))+constant_A*func1_diff(sampling_x(i));
	g2_fun(i) = func2(sampling_x(i))+constant_A*func2_diff(sampling_x(i));
	g3_fun(i) = func3(sampling_x(i))+constant_A*func3_diff(sampling_x(i));
end
%% 通过构造的反解
func1_Inverse = [f1_init];
func2_Inverse = [f2_init];
func3_Inverse = [f3_init];
for k=1:samp_num-1
	func1_Inverse(k+1) = func1_Inverse(k)+(g1_fun(k)-func1_Inverse(k))*sampling_step/constant_A;
	func2_Inverse(k+1) = func2_Inverse(k)+(g2_fun(k)-func2_Inverse(k))*sampling_step/constant_A;
	func3_Inverse(k+1) = func3_Inverse(k)+(g3_fun(k)-func3_Inverse(k))*sampling_step/constant_A;
end
%% 画图比较
raw_fun1_data = func1(sampling_x);   %准确数字
raw_fun2_data = func2(sampling_x);   %线性可以完全拟合
raw_fun3_data = func3(sampling_x);
error1 = func1_Inverse-raw_fun1_data;
error2 = func2_Inverse-raw_fun2_data;
error3 = func3_Inverse-raw_fun3_data;
%误差图
figure;
plot(sampling_x,error1,'r');
figure;
plot(sampling_x,error2,'r');
figure;
plot(sampling_x,error3,'r');
%对比图
figure;
plot(sampling_x,raw_fun1_data,'r');
hold on;
plot(sampling_x,func1_Inverse,'b');
figure;
plot(sampling_x,raw_fun2_data,'r.');
hold on;
plot(sampling_x,func2_Inverse,'b.');
figure;
plot(sampling_x,raw_fun3_data,'r.');
hold on;
plot(sampling_x,func3_Inverse,'b.');
