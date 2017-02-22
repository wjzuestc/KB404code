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
fdt=Kr*5000;
fdtt=Kr*10000;
%%
St=exp(cj*pi*Kr*t.^2);
St_OnePhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdc*t);
% St_ThreePhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdt*t.^3);
% St_fourPhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdtt*t.^4);
% St_oneAndthreePhase=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdc*t).*exp(cj*2*pi*fdtt*t.^3);
% St_final=exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdc*t).*exp(cj*2*pi*fdt*t.^3).*exp(cj*2*pi*fdtt*t.^4);
%%
figure;plot(real(St));title('LFM');
figure;plot(real(St_OnePhase));title('加一阶相');
% figure;plot(real(St_ThreePhase));title('加三阶相');
% figure;plot(real(St_fourPhase));title('加四阶相');
% figure;plot(real(St_oneAndthreePhase));title('加一阶和三阶阶相');
% figure;plot(abs(ftx(St)));title('LFM频谱');
% figure;plot(abs(ftx(St_OnePhase)));title('LFM加一阶相的频谱');
% figure;plot(abs(ftx(St_ThreePhase)));title('LFM加三阶相的频谱');
% figure;plot(real(St_final));title('LFMfinal');
% figure;plot(abs(ftx(St_final)));title('LFMfinal');
%%
% s1 = 10*St+St_OnePhase;
% s2=10*exp(cj*pi*Kr*t.^2)+exp(cj*pi*2*Kr*t.^2).*exp(cj*2*pi*0.95*Kr*t);
% s3=exp(cj*pi*Kr*t.^2)+10*exp(cj*pi*2*Kr*t.^2).*exp(cj*2*pi*0.95*Kr*t);
% figure;plot(real(s1));
% figure;plot(real(s2));
% figure;plot(real(s3));
%%
% s4=10*exp(cj*pi*Kr*t.^2).*exp(cj*2*pi*fdc*t)+exp(cj*pi*2*Kr*t.^2);
% figure;plot(real(s4));
St_final_quFdc = St_OnePhase.*exp(-cj*2*pi*fdc*t);
%%
% h = exp(-cj*pi*Kr*t.^2);
% s5 = St_final_quFdc.*h;
% figure;plot(real(s5));
% figure;plot(abs(ftx(s5)))
% [val,index] = max(abs(ftx(s5)))
%%
figure;plot(real(St_final_quFdc))
real_fdc_data = (real(St_final_quFdc)+1);
figure;plot(real_fdc_data);
%%
stepNum = size(real_fdc_data,2);
chafen = real_fdc_data(2)-real_fdc_data(1);
for i=3:stepNum
    chafenNum = real_fdc_data(i)-real_fdc_data(i-1);
    chafen = [chafen,chafenNum];
end
figure;plot(abs(chafen));

%% fdc
D=1000;
delta=fdc/D;
fdc_val = [fdc-delta*(D/2-1):delta:fdc+delta*D/2];
NumSum = 0;
% fdc=0;
%% 去模糊
for k=1:size(fdc_val,2)
St_final_quFdc = St_OnePhase.*exp(-cj*2*pi*fdc_val(k)*t);
real_fdc_data = (real(St_final_quFdc)+4)*2000;
stepNum = size(real_fdc_data,2);
%% 求一阶差分
chafen = real_fdc_data(2)-real_fdc_data(1);
for i=3:stepNum
    chafenNum = real_fdc_data(i)-real_fdc_data(i-1);
    chafen = [chafen,chafenNum];
end
% figure;plot(chafen);
%%
mid = 300;
Threshold = 50;
chafen2 = abs(chafen);
NumSum = [NumSum,sum(chafen2(mid-Threshold:mid+Threshold))];
end
%%
num_index = NumSum(2:end);
figure;
plot(fdc_val,num_index)
%%
[x,index] = min(NumSum(2:end));
fdc_estimate = fdc_val(index)
% figure;
% plot(abs(chafen));


