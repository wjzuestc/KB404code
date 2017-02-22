function Srancom_jiaozheng = Fre_domain_RCMC( varargin )
%  Ƶ�����RCMC
%  ���룺parameter--��freRCMC_Parameter��������Ĳ����б�
%        dataRR11--����ѹ����Ļز����� 
%        freRCMC_Parameter = [AziNum,fr,RanNum,fs,SapTime];
% �����Ƶ��RCMC�������
%% ��������
parameter=varargin{1};   %�����б�
Srancom=varargin{2};    %�ز�����
c=3e8;
cj=sqrt(-1);
AziNum=parameter(1);
RanNum=parameter(2);
fs=parameter(3);
SapTime=parameter(4);
PRF=parameter(5);
%%
% RanNum=3000;
%%%%%%%�����߶�У��%%%%%%%
delta_R = SapTime*c;
K = tan(3.0*pi/180)*delta_R/(PRF/AziNum);    %0.66
K1=tan(0.0033*pi/180)*delta_R/(PRF/AziNum);
% K1=tan(0*pi/180)*delta_R/(PRF/AziNum);
t0 =PRF/AziNum*K*[-AziNum/2:AziNum/2-1]/c+K1/c*[PRF/AziNum*[-AziNum/2:AziNum/2-1]].^2;
Srancom_jiaozheng=ifty((Srancom).*exp(cj*pi*2*t0.'*(fs/RanNum*[-RanNum/2:RanNum/2-1])));
end

