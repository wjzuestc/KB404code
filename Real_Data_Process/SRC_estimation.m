function [ Srancom1 ] = SRC_estimation( varargin )
%  ����SRC
%  ���룺parameter--��SRC_estimation_Parameter��������Ĳ����б�
%        dataRR11--����ѹ����Ļز����� 
%        SRC_estimation_Parameter=[BandWidth,Tpulse,AziNum];
% �����SRC����������
%% ��������
parameter=varargin{1};   %�����б�
Srancom=varargin{2};        %�ز�����
fr=varargin{3};
BandWidth=parameter(1);  %����
Tpulse=parameter(2);          %ʱ��
AziNum=parameter(3);
Kr=BandWidth/Tpulse;
%%
Filter=exp(j*pi*fr.^2/Kr/6);                  % Ƶ�����˲���
Srancom1=Srancom.*(ones(AziNum,1)*Filter);
end

