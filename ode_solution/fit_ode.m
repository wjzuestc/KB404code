% @description 用函数拟合的方式求解ode（polyfit()和ode()函数）
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
sampling_step = 0.01;         %采样间隔
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
%% 利用多项式拟合g(x)
g1_fit = polyfit(sampling_x,g1_fun,5);
g2_fit = polyfit(sampling_x,g2_fun,3);
g3_fit = polyfit(sampling_x,g3_fun,3);
g1_func_fit = @(x)g1_fit(1)*x.*x.*x.*x.*x+g1_fit(2)*x.*x.*x.*x+g1_fit(3)*x.*x.*x+g1_fit(4)*x.*x+g1_fit(5)*x+g1_fit(6);
g2_func_fit = @(x)g2_fit(1)*x.*x.*x+g2_fit(2)*x.*x+g2_fit(3)*x+g2_fit(4);
g3_func_fit = @(x)g3_fit(1)*x.*x.*x+g3_fit(2)*x.*x+g3_fit(3)*x+g3_fit(4);
%% ode求解
% odefun1 = @(x,f1)(g1_fit(1)*x.*x.*x+g1_fit(2)*x.*x+g1_fit(3)*x+g1_fit(4)-f1)/constant_A;
odefun1 = @(x,f1)(g1_fit(1)*x.*x.*x.*x.*x+g1_fit(2)*x.*x.*x.*x+g1_fit(3)*x.*x.*x+g1_fit(4)*x.*x+g1_fit(5)*x+g1_fit(6))/constant_A;
odefun2 = @(x,f2)(g2_fit(1)*x.*x.*x+g2_fit(2)*x.*x+g2_fit(3)*x+g2_fit(4)-f2)/constant_A;
odefun3 = @(x,f3)(g3_fit(1)*x.*x.*x+g3_fit(2)*x.*x+g3_fit(3)*x+g3_fit(4)-f3)/constant_A;
f1_init = func1(0.01);           %设置初始点
f2_init = func2(0.01);
f3_init = func3(0.01);
[x,f1] = ode45(odefun1,[0.01:0.01: 10],f1_init);
[x,f2] = ode45(odefun2,[0.01:0.01: 10],f2_init);
[x,f3] = ode45(odefun3,[0.01:0.01: 10],f3_init);
%% 画图比较
sampling_x = [0.01:0.01: 10];
raw_fun1_data = func1(sampling_x);      %准确数字
raw_fun2_data = func2(sampling_x);      %线性可以完全拟合
raw_fun3_data = func3(sampling_x);
error1 = f1'-raw_fun1_data;
error2 = f2'-raw_fun2_data;
error3 = f3'-raw_fun3_data;
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
plot(x,f1,'b');
figure;
plot(sampling_x,raw_fun2_data,'r.');
hold on;
plot(x,f2,'b.');
figure;
plot(sampling_x,raw_fun3_data,'r.');
hold on;
plot(x,f3,'b.');
