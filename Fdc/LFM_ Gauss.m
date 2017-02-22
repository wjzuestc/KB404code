clc;clear;close all;
%% ����
cj=sqrt(-1);
c=3e8;
%% ϵͳ����
Tp=10e-6;   %ʱ��                               
BandWidth=30e6;    %����                         
Kr=BandWidth/Tp;   %��Ƶб��                              
Fs=2*BandWidth;    %������
Ts=1/Fs;           %����ʱ��            
N=Tp/Ts;           %��������
t=linspace(-Tp/2,Tp/2,N);
%Kr 3e12
fdc=0.9*Kr;
fdt=Kr*5000;
fdtt=Kr*10000;
%%
Gauss=wgn(1,600,10,'dBm','complex');
figure;plot(real(Gauss));title('Gauss�������ź�')
%%
St=exp(cj*pi*Kr*t.^2);
St_OnePhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdc*t);
St_ThreePhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdt*t.^3);
St_fourPhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdtt*t.^4);
St_oneAndthreePhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdc*t).*exp(cj*2*pi*fdtt*t.^3);
St_final=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdc*t).*exp(cj*2*pi*fdt*t.^3).*exp(cj*2*pi*fdtt*t.^4);
St_Gauss=St+Gauss;
St_OnePhase_Gauss=St_OnePhase+Gauss;
% St_Gauss=awgn(St,'SNR','measured')
%%
figure;plot(real(St_Gauss));title('LFM����Gauss�ź�')
figure;plot(real(St_OnePhase_Gauss));title('LFM�źż���fdc��GaussӰ��')
%%
% figure;plot(real(St));title('LFM');
% figure;plot(real(St_OnePhase));title('��һ����');
% figure;plot(real(St_ThreePhase));title('��������');
% figure;plot(real(St_fourPhase));title('���Ľ���');
% figure;plot(real(St_oneAndthreePhase));title('��һ�׺����׽���');
% figure;plot(abs(ftx(St)));title('LFMƵ��');
% figure;plot(abs(ftx(St_OnePhase)));title('LFM��һ�����Ƶ��');
% figure;plot(abs(ftx(St_ThreePhase)));title('LFM���������Ƶ��');
% figure;plot(real(St_final));title('LFMfinal');
% figure;plot(abs(ftx(St_final)));title('LFMfinal');
%%
St_OnePhase_Gauss_quFdc=St_OnePhase_Gauss.*exp(-cj*2*pi*fdc*t) 
figure;plot(real(St_OnePhase_Gauss_quFdc));
stepNum = size(St_OnePhase_Gauss_quFdc,2);
chafen = St_OnePhase_Gauss_quFdc(2)-St_OnePhase_Gauss_quFdc(1);
for i=3:stepNum
    chafenNum = St_OnePhase_Gauss_quFdc(i)-St_OnePhase_Gauss_quFdc(i-1);
    chafen = [chafen,chafenNum];
end
figure;plot(abs(chafen));title('ȥģ�����һ��΢��')
