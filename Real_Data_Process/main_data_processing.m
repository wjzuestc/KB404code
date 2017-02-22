%% ʵ�����ݴ�����ֹƽ̨��
clc;
close all;
clear all;
%% ��Ҫ����
Nrange=30000;           %�������������
Nblock=10000;
fs=1e9;                 %�����������
Tsamplemin=523/3e8;     %������ʼʱ��
Tpulse=2e-5;            %����
BandWidth=300e6;        %����
PRF=2000;               %�������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���ù�
Nazimuth=4000;          %��λ�������
Nchou=1;               
Nmin=1;                 %�������ݾ��������
Nmax=Nrange;            %�������ݾ������յ�
Nframe=6;
Nhangzi=94;
fs=1e9/Nchou;
Rbin=3e8/fs;
Rmin=Tsamplemin*3e8-Tpulse*3e8;                              %
Rmax=Rmin+(Nmax-Nmin)*Rbin;
fc = 750e6;                %��Ƶ
Nrange=Nrange+Nhangzi+Nframe;   %6֡ͷ+32����+�ز�
Rrange=linspace(Rmin,Rmax,Nmax-Nmin+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filePath=pwd;                  %�ز��ļ�·��
first=1;
% last=1;
last=15;
AziNum=Nazimuth*(last-first+1);
commonParameter=[0;BandWidth;Tpulse;1/PRF;Nazimuth;fc;fs;Nblock]; %[0������ʱ���������ʱ�䣬��λ������㣬��Ƶ�����������Ƶ�ʣ�]
%% *******************������************************
dataRR=[];              %�ز�����
for dijisao = first:last
    [dataR,echoData]=ReaddataFMF(commonParameter,0,dijisao,Nrange,1,filePath,Nmin,Nmax,1);%���ز�    
    dataRR=[dataRR,dataR];   
end
%% *******************����������ѹ��************************
% ���Ӵ�
range_compression_Parameter=[BandWidth;Tpulse;AziNum;fs];
[dataRRR]=range_compression(range_compression_Parameter,dataRR); 
%�Ӵ�
% Nrange=30000;
% [dataRRR]=range_compression_windows(range_compression_Parameter,dataRR,Nrange); 
% figure
% imagesc((abs((dataRRR.'))))
%%
%ȡ����
a=23300;
RanNum=2000;
dataRR11=dataRRR(a-RanNum/2+1:a+RanNum/2,:);  
figure
imagesc((abs((dataRR11.'))))
%%
SapTime = 1/fs; 
delta_R = SapTime*3e8;
Nrange=2000;
PRI=1/PRF;
Ta = PRI*[-AziNum/2:AziNum/2-1];                       %��ʱ��
fc=9.6e9;
c=3e8;
lambda=c/fc;
cj=sqrt(-1);
RT=430;RR=475;
R0=RT+RR;
V=30;
Ts=R0*(9/180*pi)/V;
phiR=25/180*pi;
x=200/tan(phiR)-280;
phiT=asin(x/RT);
fdc=V*sin(phiT)/lambda+V*cos(phiR)/lambda;
fdr=-V^2*cos(phiT)/lambda/RT-V^2*sin(phiR)/lambda/RR;
fr=fs/RanNum*[-RanNum/2:RanNum/2-1];                  %������ο�Ƶ��
fa=PRF/AziNum*[-AziNum/2:AziNum/2-1]+fdc;             %������ο�Ƶ��
%% *******************RCMC************************
%ʱ��RCMC
timeRCMC_Parameter = [fdc,lambda,V,fdr];
S_shiyuRCMC=time_domain_RCMC(timeRCMC_Parameter,dataRR11,fr,Ta);
figure;
imagesc(abs(S_shiyuRCMC))
Sran=iftx(S_shiyuRCMC);

%% PRF��бУ��
dataRR1=dataRR11.*(ones(RanNum,1)*exp(cj*pi*(fdc-2200)*Ta));
S=dataRR1.';
PRF_Correction_Parameter=[AziNum,RanNum];
S_PRF_Correction = PRF_Tilt_correction(PRF_Correction_Parameter,S,fr);
figure;colormap(gray);imagesc(abs(S_PRF_Correction))
%% Ƶ��RCMC
Srancom=ftx(fty(S_PRF_Correction));
freRCMC_Parameter = [AziNum,RanNum,fs,SapTime,PRF];
S_pingyuRCMC=Fre_domain_RCMC(freRCMC_Parameter,Srancom);
SranPingyu=iftx(S_pingyuRCMC);
figure
imagesc(abs(SranPingyu))
%% *********************SRC**************************
%����SRC
% SRC_estimation_Parameter=[BandWidth,Tpulse,AziNum];
% Srancom1 = SRC_estimation(SRC_estimation_Parameter,Sran,fr);
% figure;imagesc(abs(ifty(Srancom1)))
% figure;imagesc(abs(ifty(Sran)))
%����SRC
SRC_lilun_Parameter=[a,Nrange,fc,SapTime];
Sran = lilun_SRC(SRC_lilun_Parameter,SranPingyu,Rrange,fa,fr);
figure;imagesc(abs(ifty(Sran)))
%% *******************��λ��ѹ��************************
Imageraw = azimuth_compression(fdr,Ta,Ts,Nrange,Sran);




