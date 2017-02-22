function [dataout] = range_compression( varargin )
%  �Իز����о���������ѹ�������Ӵ���
%  ���룺parameter--��Readparameter��������Ĳ����б�
%        data--�ز����� 
%        range_compression_Parameter=[BandWidth;Tpulse;AziNum;fs];
% ���������������ѹ����Ļز�����
%% ��������
parameter=varargin{1};   %�����б�
data=varargin{2};        %�ز�����
Bandwidth=parameter(1);  %����
T=parameter(2);          %ʱ��
Nazimuth=parameter(3);
samplerate_full = parameter(4);

sample_beishu = 1;   %������
cj=sqrt(-1);
fs=samplerate_full/sample_beishu; % ������
Kr=Bandwidth/T;                   %��Ƶб��
[Mdata,Ndata]=size(data); 
Nrange=Mdata;
%% ����������ѹ�����Ӵ�
t1=linspace(-T/2,T/2,T*fs);
Filter=exp(-cj*pi*Kr*t1.^2);                  % ʱ�����˲���
TTT=Nrange-length(Filter);
Filter=[Filter zeros(1,TTT)];
data = fft(data).*fft(Filter.'*ones(1,Nazimuth));     % Ƶ���˲�
data = ifft(data);

dataout=data;
end

