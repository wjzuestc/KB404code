clc;clear all;close all;
%% 参数
cj=sqrt(-1);
c=3e8;
%% % (I) parameters' definition
res_a=2;				    % 所需的方位分辨率
res_r=2;				    % 所需的距离分辨率
k_a=1.2;				    % 方位因子
k_r=1.2;				    % range factor

Ra=5000.;				    % 雷达工作距离
va=70.; 				    % radar/platform forward velocity速度
Tp=1.e-6;					% transmitted pulse width 脉宽     
fc=3.e+9;	 				% carrier frequency 频率
FsFactor = 1.2;
theta=90*pi/180;			% squint angle斜视角  

lamda=c/fc;							    % 波长
Br=k_r*c/2./res_r;					    % required transmitted bandwidth传输带宽
Fs=Br*FsFactor;						    % A/D sampling rate采样频率
Kr=Br/Tp;					  		    %调频斜率	  

La=Ra*k_a*lamda/2/res_a/sin(theta);     % 合成孔径长度
Ta=La/va;							    % 合成孔径时间
fdc=2*va*cos(theta)/lamda;              % doppler centriod
fdr=-2*(va*sin(theta)).^2/lamda/Ra;	    % doppler rate
Bd=abs(fdr)*Ta;						    % doppler bandwidth
prf=round(Bd*2);					    % PRF脉冲重复频率	
%% =====================================================================                                   
Nr=Tp*Fs;
Ni=1:Nr;
tr=(Ni-Nr/2)*Tp/Nr;
St=exp(j*pi*Kr*(tr).^2).*exp(j*2*pi*2e13*tr).*exp(j*2*pi*2e13*(tr).^3);                   %生成线性调频信号
figure;plot(real(St))
figure;plot(abs(ftx(St)))