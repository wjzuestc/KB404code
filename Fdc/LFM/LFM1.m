clc;clear all;close all;
%% 参数
cj=sqrt(-1);
c=3e8;
%%
%系统参数
Location_T=[-1000,0,1000];%发射机平台的初始位置，向量第一项对应距离，第二项对应方位，第三项对应高度
Location_R=[0,-500,1000];%接收平台初始时刻的位置
Vt=[0 50 0];%接收站速度
Vr=[0 50 0];%发射站速度
%% 信号参数
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
%% calculate value
TargetCenter=[0 0 0];%目标中心
Rsq_R=sqrt((TargetCenter-Location_R)*(TargetCenter-Location_R).');
Rsq_T=sqrt((TargetCenter-Location_T)*(TargetCenter-Location_T).');
R_center=Rsq_T+Rsq_R;

Tr=SapTime*(-RanNum/2:RanNum/2-1);  %距离向时间
Ta=PRI*(-AziNum/2:AziNum/2-1);      %方位向时间
T_path=ones(AziNum,1)*Location_T+Ta.'*Vt;
R_path=ones(AziNum,1)*Location_R+Ta.'*Vr;
%%
fdc = 789;
Targets=[0 0 0];
Targets=Targets+ones(size(Targets,1),1)*TargetCenter;
RawData=zeros(AziNum,RanNum);
for i=1:size(Targets,1)
    AziNum_Rr=sqrt(sum((R_path-ones(AziNum,1)*Targets(i,:)).^2,2));
    AziNum_Rt=sqrt(sum((T_path-ones(AziNum,1)*Targets(i,:)).^2,2));
    R_target=AziNum_Rt+AziNum_Rr;
    td=ones(AziNum,1)*Tr-(R_target-R_center)/c*ones(1,RanNum); %t-时延   
%     RawData=RawData+exp(cj*2*pi*2.0e13*td).*exp(-cj*2*pi*fc*R_target/c*ones(1,RanNum)+cj*pi*Kr*td.^2).*...
%         (abs(td)<Tp/2).*(abs(Targets(i,2)-R_path(:,2)-(TargetCenter(2)-Location_R(2)))*ones(1,RanNum)<Ts*Vr(2)/2);   
    RawData=RawData+exp(cj*2*pi*2.0e13*td).*exp(cj*pi*Kr*td.^2).*(abs(td)<Tp/2);   
end
figure;
imagesc(real(RawData(1,:)));
figure('Name','回波');
imagesc(abs(ftx(RawData(1,:)))); %频谱图