function  S1  = PRF_Tilt_correction( varargin )
%  PRF��бУ��
%  ���룺parameter--��freRCMC_Parameter��������Ĳ����б�
%        dataRR11--����ѹ����Ļز����� 
%        freRCMC_Parameter = [AziNum,fr,RanNum,fs,SapTime];
% �����PRFУ���������
%% ��������
parameter=varargin{1};   %�����б�
S=varargin{2};    %�ز�����
fr=varargin{3};
c=3e8;
cj=sqrt(-1);
AziNum=parameter(1);
RanNum=parameter(2);
%% PRF��бУ��
Fr=ones(AziNum,1)*fr;
Rerror = linspace(0,123.84,AziNum).'*ones(1,RanNum);        %%ʵ�����������������������ĵľ�����ʷ��֮��     
S1 =ifty(fty(S).*exp(cj*2*pi/c*Rerror.*Fr));        %%�������ز���������ز��洢λ�ô�λ��Խ���뵥Ԫ��
% S2=S1(:,1001:4000);
end

