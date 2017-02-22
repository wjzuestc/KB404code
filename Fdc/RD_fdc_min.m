clc;clear all;close all;
%%
%参数
cj=sqrt(-1);
c=3e8;
%%
%系统参数
Location_T=[-1000,0,1000];%发射机平台的初始位置，向量第一项对应距离，第二项对应方位，第三项对应高度
Location_R=[0,-500,1000];%接收平台初始时刻的位置
Vt=[0 50 0];%接收站速度
Vr=[0 50 0];%发射站速度
%%
%信号参数
Tp=8E-6;%发射信号的时宽
BandWidth=200E6;%发射信号的带宽
Kr = BandWidth/Tp;%发射信号的调频斜率
SapRate = 1.5*BandWidth;%距离向采样频率
SapTime = 1/SapRate;%距离向采样时间间隔
fc = 9.65E+9;%载波频率
lambda = c/fc;%载波波长
PRF = 390;%脉冲重复频率
PRI = 1/PRF;%脉冲重复时间
Ts=1.8;%方位向照射时间
AziNum = 2048*2;%方位向采样点数
RanNum = 2048*2;%距离向采样点数
%%
%calculate value
TargetCenter=[0 0 0];%目标中心
Rsq_R=sqrt((TargetCenter-Location_R)*(TargetCenter-Location_R).');
Rsq_T=sqrt((TargetCenter-Location_T)*(TargetCenter-Location_T).');
R_center=Rsq_T+Rsq_R;

Tr=SapTime*(-RanNum/2:RanNum/2-1);  %距离向时间
Ta=PRI*(-AziNum/2:AziNum/2-1);     %方位向时间
T_path=ones(AziNum,1)*Location_T+Ta.'*Vt;
R_path=ones(AziNum,1)*Location_R+Ta.'*Vr;
%%
%create echo
% Targets=[0 0 0;0 255 0;0 -255 0;255 0 0;-255 0 0;255 255 0;-255 255 0;-255 -255 0;255 -255 0];
Targets=[0 0 0];
Targets=Targets+ones(size(Targets,1),1)*TargetCenter;
RawData=zeros(AziNum,RanNum);
for i=1:size(Targets,1)
    i
    AziNum_Rr=sqrt(sum((R_path-ones(AziNum,1)*Targets(i,:)).^2,2));
    AziNum_Rt=sqrt(sum((T_path-ones(AziNum,1)*Targets(i,:)).^2,2));
    R_target=AziNum_Rt+AziNum_Rr;
    td=ones(AziNum,1)*Tr-(R_target-R_center)/c*ones(1,RanNum); %t-时延
    
    RawData=RawData+exp(-cj*2*pi*fc*R_target/c*ones(1,RanNum)+cj*pi*Kr*td.^2).*...
        (abs(td)<Tp/2).*(abs(Targets(i,2)-R_path(:,2)-(TargetCenter(2)-Location_R(2)))*ones(1,RanNum)<Ts*Vr(2)/2);   
end
figure;
imagesc(real(RawData));colormap(gray);
figure('Name','回波');
imagesc(abs(fty(ftx(RawData)))); %频谱图
%%
%距离向脉压
hr=fftshift((abs(Tr)<Tp/2).*(0.54-0.46*cos(2*pi*(1:RanNum)/RanNum)).*...
        exp(cj*pi*Kr*(Tr).^2));
Srancom=fft(RawData.').*conj(fft(hr.')*ones(1,AziNum));
Srancom=ifft(Srancom).';
figure;imagesc(abs(Srancom));colormap(gray);
%%
thetaT=asin(Vr*(TargetCenter-Location_T).'/sqrt(Vt*Vt.')/...
    sqrt((TargetCenter-Location_T)*(TargetCenter-Location_T).'));
thetaR=asin(Vr*(TargetCenter-Location_R).'/sqrt(Vr*Vr.')/...
    sqrt((TargetCenter-Location_R)*(TargetCenter-Location_R).'));
Referent_Rsq_T=norm(TargetCenter-Location_T);
Referent_Rsq_R=norm(TargetCenter-Location_R);

k0=Referent_Rsq_T+Referent_Rsq_R;
k1=-Vt(2)*sin(thetaT)-Vr(2)*sin(thetaR);
k2=Vt(2)^2*cos(thetaT)^2/Referent_Rsq_T+Vr(2)^2*cos(thetaR)^2/Referent_Rsq_R;
k3=3*Vt(2)^3*cos(thetaT)^2*sin(thetaT)/Referent_Rsq_T^2+3*Vr(2)^3*cos(thetaR)^2*sin(thetaR)/Referent_Rsq_R^2;
k4=3*Vt(2)^4*cos(thetaT)^2*(5*sin(thetaT)^2-1)/Referent_Rsq_T^3+3*Vr(2)^4*cos(thetaR)^2*(5*sin(thetaR)^2-1)/Referent_Rsq_R^3;

theta_sqx=asin(k3/sqrt(5*k3^2-3*k2*k4));
Vx=k2*k3*(5*sin(theta_sqx)^2-1)/(2*k4*cos(theta_sqx)^2*sin(theta_sqx));
Rx=3*k2*Vx*sin(theta_sqx)/k3;
a0=k0/2-Rx;
a1=k1/2-Vx*sin(theta_sqx);

referent_fdc=sqrt(Vr*Vr.')*sin(thetaR)/lambda+...
    sqrt(Vt*Vt.')*sin(thetaT)/lambda;

ft=SapRate/RanNum*[-RanNum/2:RanNum/2-1];%距离向参考频率
fs=PRF/AziNum*[-AziNum/2:AziNum/2-1];%方位向参考频率fdc
fs_field=round(referent_fdc/PRF);
fdc_index=round((referent_fdc-PRF*fs_field)/(PRF/AziNum));
if fdc_index>=0
    fs(1,fdc_index:end)=fs(1,fdc_index:end)+fs_field*PRF;
    fs(1,1:fdc_index-1)=fs(1,1:fdc_index-1)+(fs_field+1)*PRF;
else
    fs(1,1:AziNum+fdc_index)=fs(1,1:AziNum+fdc_index)+fs_field*PRF;
    fs(1,AziNum+fdc_index+1:AziNum)=fs(1,AziNum+fdc_index+1:AziNum)+(fs_field-1)*PRF;
end
% figure;plot(fs);
a1=0;
[ft_m,fs_m]=meshgrid(ft,fs);
formulaT=sqrt((fc/c)^2-((fs_m*c+2*a1*fc)/(2*Vt(2)*c)).^2);
formulaR=sqrt((fc/c)^2-((fs_m*c+2*a1*fc)/(2*Vr(2)*c)).^2);
formula=fc/c^2-a1*(c*fs_m+2*a1*fc)/2/Vt(2)^2/c^2;
%%
%成像算法
SS=fftshift(fft(fft(Srancom).').');
%SRC
theta_src=-2*pi*Rsq_T*cos(thetaT)*((1/c^2-a1^2/Vt(2)^2/c^2)*formulaT.^2-formula.^2)./formulaT.^3.*ft_m.^2+...
    -2*pi*Rsq_R*cos(thetaR)*((1/c^2-a1^2/Vr(2)^2/c^2)*formulaR.^2-formula.^2)./formulaR.^3.*ft_m.^2;
%RCM
theta_rcm=-2*pi*Rsq_T*cos(thetaT)*formula./formulaT-...
    2*pi*Rsq_R*cos(thetaR)*formula./formulaR;
theta_rcm=(theta_rcm-theta_rcm(AziNum/2+1,RanNum/2+1)).*ft_m;
SS=SS.*exp(-cj*theta_rcm);
%二维逆傅里叶变换
sss=ifft(ifft(fftshift(SS)).').';
figure;imagesc(abs(sss));colormap(gray);
%% fdc
fdc = [600:10:850];
NumSum = 0;
% fdc=0;
%% 去模糊
for k=1:size(fdc,2)
ss=(sss.'.*exp(-1i*2*pi*fdc(k)*ones(RanNum,1)*Ta)).';
% figure;imagesc(real(ss));colormap(gray);
% figure('Name','去模糊的回波');
% imagesc(abs(fty(ftx(ss))));
%% 取一个距离门
[val,RanSelect] = max(max(real(ss)));
fdc_process_data=ss(1699:2399,RanSelect);
% figure;plot(imag(fdc_process_data))
% figure;plot(real(fdc_process_data))
%%
real_fdc_data = real(fdc_process_data)+2000;
% figure;plot(real_fdc_data);
stepNum = size(real_fdc_data,1);
%% 求一阶差分
chafen = real_fdc_data(2)-real_fdc_data(1);
for i=3:stepNum
    chafenNum = real_fdc_data(i)-real_fdc_data(i-1);
    chafen = [chafen,chafenNum];
end
% figure;plot(chafen);
%%
mid = 350;
chafen2 = abs(chafen);
NumSum = [NumSum,sum(chafen2(mid-10:mid+10))];
end
%%
figure;plot(fdc,NumSum(2:end));xlabel('fdc值');ylabel('Z(fdc)')
%%
[x,index] = min(NumSum(2:end));
fdc_estimate = fdc(index)




