clc;clear;close all;
%% 参数
cj=sqrt(-1);
c=3e8;
%% 系统参数
Tp=10e-6;   %时宽                               
BandWidth=30e6;    %带宽                         
Kr=BandWidth/Tp;   %调频斜率                              
Fs=2*BandWidth;    %采样率
Ts=1/Fs;           %采样时间            
N=Tp/Ts;           %采样点数
t=linspace(-Tp/2,Tp/2,N);
%Kr 3e12
fdc=0.9*Kr;
fdt=Kr*1e9;
fdtt=Kr*10000;
%%
St=exp(cj*pi*Kr*t.^2);
St_OnePhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdc*t);
St_ThreePhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdt*t.^3);
St_fourPhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdtt*t.^4);
St_oneAndthreePhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdc*t).*exp(cj*2*pi*fdtt*t.^3);
% St_final=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdc*t).*exp(cj*2*pi*fdt*t.^3).*exp(cj*2*pi*fdtt*t.^4);
%%
figure;plot(real(St));title('LFM');
figure;plot(real(St_OnePhase));title('加一阶相');
figure;plot(real(St_ThreePhase));title('加三阶相');
figure;plot(real(St_fourPhase));title('加四阶相');
figure;plot(real(St_oneAndthreePhase));title('加一阶和三阶阶相');
% figure;plot(abs(ftx(St)));title('LFM频谱');
% figure;plot(abs(ftx(St_OnePhase)));title('LFM加一阶相的频谱');
% figure;plot(abs(ftx(St_ThreePhase)));title('LFM加三阶相的频谱');
% figure;plot(real(St_final));title('LFMfinal');
% figure;plot(abs(ftx(St_final)));title('LFMfinal');
