clc;clear;close all;
%% 参数
cj=sqrt(-1);
c=3e8;
%% 系统参数
Tp=1;   %时宽                               
BandWidth=1000;    %带宽                         
Kr=BandWidth/Tp;   %调频斜率                              
Fs=2*BandWidth;    %采样率
Ts=1/Fs;           %采样时间            
N=Tp/Ts;           %采样点数
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
figure;plot(real(St_OnePhase));title('加一阶相');
figure;plot(real(St_ThreePhase));title('加三阶相');
figure;plot(real(St_fourPhase));title('加四阶相');
% figure;plot(real(St_oneAndthreePhase));title('加一阶和三阶阶相');
% figure;plot(abs(ftx(St)));title('LFM频谱');
% figure;plot(abs(ftx(St_OnePhase)));title('LFM加一阶相的频谱');
% figure;plot(abs(ftx(St_ThreePhase)));title('LFM加三阶相的频谱');
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
% figure;plot(chafen);title('一阶微分')

%% fdc
fdc_val = [500000:10000:1500000];
NumSum = 0;
% fdc=0;
%% 去模糊
for k=1:size(fdc_val,2)
St_final_quFdc = St_final.*exp(-cj*2*pi*fdc_val(k)*t);
real_fdc_data = (real(St_final_quFdc)+1)*2000;
stepNum = size(real_fdc_data,2);
%% 求一阶差分
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
figure;plot(fdc_val,num_index);title('Z(fdc)');xlabel('fdc值');ylabel('Z(fdc)值')
%%
[x,index] = min(NumSum(2:end));
fdc_estimate = fdc_val(index)
% percent = (fdc_estimate-fdc)/fdc



