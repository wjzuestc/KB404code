clc;clear all;close all;
%% ����
cj=sqrt(-1);
c=3e8;
%%
%ϵͳ����
Location_T=[-1000,0,1000];%�����ƽ̨�ĳ�ʼλ�ã�������һ���Ӧ���룬�ڶ����Ӧ��λ���������Ӧ�߶�
Location_R=[0,-500,1000];%����ƽ̨��ʼʱ�̵�λ��
Vt=[0 50 0];%����վ�ٶ�
Vr=[0 50 0];%����վ�ٶ�
%% �źŲ���
Tp=8E-6;%�����źŵ�ʱ��
BandWidth=200E6;%�����źŵĴ���
Kr = BandWidth/Tp;%�����źŵĵ�Ƶб��
SapRate = 1.5*BandWidth;%���������Ƶ��
SapTime = 1/SapRate;%���������ʱ����
fc = 9.65E+9;%�ز�Ƶ��
lambda = c/fc;%�ز�����
PRF = 390;%�����ظ�Ƶ��
PRI = 1/PRF;%�����ظ�ʱ��
Ts=1.8;%��λ������ʱ��
AziNum = 2048*2;%��λ���������
RanNum = 2048*2;%�������������
%% calculate value
TargetCenter=[0 0 0];%Ŀ������
Rsq_R=sqrt((TargetCenter-Location_R)*(TargetCenter-Location_R).');
Rsq_T=sqrt((TargetCenter-Location_T)*(TargetCenter-Location_T).');
R_center=Rsq_T+Rsq_R;

Tr=SapTime*(-RanNum/2:RanNum/2-1);  %������ʱ��
Ta=PRI*(-AziNum/2:AziNum/2-1);      %��λ��ʱ��
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
    td=ones(AziNum,1)*Tr-(R_target-R_center)/c*ones(1,RanNum); %t-ʱ��   
%     RawData=RawData+exp(cj*2*pi*2.0e13*td).*exp(-cj*2*pi*fc*R_target/c*ones(1,RanNum)+cj*pi*Kr*td.^2).*...
%         (abs(td)<Tp/2).*(abs(Targets(i,2)-R_path(:,2)-(TargetCenter(2)-Location_R(2)))*ones(1,RanNum)<Ts*Vr(2)/2);   
    RawData=RawData+exp(cj*2*pi*2.0e13*td).*exp(cj*pi*Kr*td.^2).*(abs(td)<Tp/2);   
end
figure;
imagesc(real(RawData(1,:)));
figure('Name','�ز�');
imagesc(abs(ftx(RawData(1,:)))); %Ƶ��ͼ