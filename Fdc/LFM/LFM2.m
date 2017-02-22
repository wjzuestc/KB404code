clc;clear all;close all;
%% ����
cj=sqrt(-1);
c=3e8;
%% % (I) parameters' definition
res_a=2;				    % ����ķ�λ�ֱ���
res_r=2;				    % ����ľ���ֱ���
k_a=1.2;				    % ��λ����
k_r=1.2;				    % range factor

Ra=5000.;				    % �״﹤������
va=70.; 				    % radar/platform forward velocity�ٶ�
Tp=1.e-6;					% transmitted pulse width ����     
fc=3.e+9;	 				% carrier frequency Ƶ��
FsFactor = 1.2;
theta=90*pi/180;			% squint angleб�ӽ�  

lamda=c/fc;							    % ����
Br=k_r*c/2./res_r;					    % required transmitted bandwidth�������
Fs=Br*FsFactor;						    % A/D sampling rate����Ƶ��
Kr=Br/Tp;					  		    %��Ƶб��	  

La=Ra*k_a*lamda/2/res_a/sin(theta);     % �ϳɿ׾�����
Ta=La/va;							    % �ϳɿ׾�ʱ��
fdc=2*va*cos(theta)/lamda;              % doppler centriod
fdr=-2*(va*sin(theta)).^2/lamda/Ra;	    % doppler rate
Bd=abs(fdr)*Ta;						    % doppler bandwidth
prf=round(Bd*2);					    % PRF�����ظ�Ƶ��	
%% =====================================================================                                   
Nr=Tp*Fs;
Ni=1:Nr;
tr=(Ni-Nr/2)*Tp/Nr;
St=exp(j*pi*Kr*(tr).^2).*exp(j*2*pi*2e13*tr).*exp(j*2*pi*2e13*(tr).^3);                   %�������Ե�Ƶ�ź�
figure;plot(real(St))
figure;plot(abs(ftx(St)))