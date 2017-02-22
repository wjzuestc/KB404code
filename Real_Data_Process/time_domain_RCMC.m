function [ Srancom_jiaozheng ] = time_domain_RCMC( varargin )
%  ʱ�����RCMC
%  ���룺parameter--��timeRCMC_Parameter��������Ĳ����б�
%        dataRR11--����ѹ����Ļز����� 
%        timeRCMC_Parameter = [fdc,lambda,V,Ta,fr,fdr];
% �����ʱ��RCMC�������
%% ��������
parameter=varargin{1};   %�����б�
dataRR11=varargin{2};        %�ز�����
fr=varargin{3};
Ta=varargin{4};
c=3e8;
cj=sqrt(-1);
fdc=parameter(1);
lambda=parameter(2);
V=parameter(3);
fdr=parameter(4);
%%
Sy=fty(dataRR11.');
K0=fdc/1.5*lambda/V;%%%%�����㶯У����Ҫ���ԵĲ���%%%%%%%%
t_ij=-V*(K0)*Ta/c;
comf=exp(cj*2*pi*t_ij.'*fr);             
Sy1=Sy.*comf;
S1=ifty(Sy1);


K1=fdr/8*lambda;%%%%�����㶯У����Ҫ���ԵĲ���%%%%%%%%
t_ij1=-(K1)*Ta.^2/c;
comf=exp(cj*2*pi*t_ij1.'*fr);             
Sy2=Sy1.*comf;
Srancom_jiaozheng=ifty(Sy2);
end

