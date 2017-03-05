clc;clear;close all;
%% ����
cj=sqrt(-1);
c=3e8;
%% ϵͳ����
Tp=1;   %ʱ��                               
BandWidth=1000;    %����                         
Kr=BandWidth/Tp;   %��Ƶб��                              
Fs=2*BandWidth;    %������
Ts=1/Fs;           %����ʱ��            
N=Tp/Ts;           %��������
t=linspace(-Tp/2,Tp/2,N);
%Kr 3e12
fdc=1000*Kr;
fdt=Kr*5000;
fdtt=Kr*10000;
%%
St=exp(cj*pi*Kr*t.^2);
St_OnePhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdc*t);
St_ThreePhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdt*t.^3);
St_fourPhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdtt*t.^4);
% St_oneAndthreePhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdc*t).*exp(cj*2*pi*fdtt*t.^3);
St_final=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdc*t).*exp(cj*2*pi*fdt*t.^3).*exp(cj*2*pi*fdtt*t.^4);
%%
figure;plot(real(St));title('LFM');
figure;plot(real(St_OnePhase));title('��һ����');
figure;plot(real(St_ThreePhase));title('��������');
figure;plot(real(St_fourPhase));title('���Ľ���');
% figure;plot(real(St_oneAndthreePhase));title('��һ�׺����׽���');
% figure;plot(abs(ftx(St)));title('LFMƵ��');
% figure;plot(abs(ftx(St_OnePhase)));title('LFM��һ�����Ƶ��');
% figure;plot(abs(ftx(St_ThreePhase)));title('LFM���������Ƶ��');
figure;plot(real(St_final));title('LFMfinal');
figure;plot(abs(ftx(St_final)));title('LFMfinal');
%%
% St_final_quFdc = St_final.*exp(-cj*2*pi*fdc*t);
% figure;plot(real(St_final_quFdc))
% %%
% real_fdc_data = (real(St_final_quFdc)+1)*2000;
% figure;plot(real_fdc_data);
% stepNum = size(real_fdc_data,2);
% chafen = real_fdc_data(2)-real_fdc_data(1);
% for i=3:stepNum
%     chafenNum = real_fdc_data(i)-real_fdc_data(i-1);
%     chafen = [chafen,chafenNum];
% end
% figure;plot(chafen);title('һ��΢��')

%% fdc
fdc_val = [500000:10000:1500000];
NumSum = 0;
% fdc=0;
%% ȥģ��
for k=1:size(fdc_val,2)
St_final_quFdc = St_final.*exp(-cj*2*pi*fdc_val(k)*t);
real_fdc_data = (real(St_final_quFdc)+1)*2000;
stepNum = size(real_fdc_data,2);
%% ��һ�ײ��
chafen = real_fdc_data(2)-real_fdc_data(1);
for i=3:stepNum
    chafenNum = real_fdc_data(i)-real_fdc_data(i-1);
    chafen = [chafen,chafenNum];
end
% figure;plot(chafen);
%%
mid = 1000;
Threshold = 2;
chafen2 = abs(chafen);
NumSum = [NumSum,sum(chafen2(mid-Threshold:mid+Threshold))];
end
%%
num_index = NumSum(2:end);
figure;plot(fdc_val,num_index);title('Z(fdc)');xlabel('fdcֵ');ylabel('Z(fdc)ֵ')
%%
[x,index] = min(NumSum(2:end));
fdc_estimate = fdc_val(index)
% percent = (fdc_estimate-fdc)/fdc



