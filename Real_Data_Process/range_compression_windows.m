function [ dataout ] = range_compression_windows( varargin )
%  �Իز����о���������ѹ��--�Ӵ�
%  ���룺parameter--��Readparameter��������Ĳ����б�
%        data--�ز����� 
%        Nrange=varargin{3};
% parameter=[fuyangjiao;Bandwidth;T;PRF;Vscan;Amin;Amax;H;Va;AngleVH;AH;AV;AR;Longitude;Latitude;dushift;Nazimuth];
% ���������������ѹ�����Ӵ�����Ļز�����
%% ��������
parameter=varargin{1};     %�����б�
dataRR=varargin{2};        %�ز�����
Nrange=varargin{3};        %�������������
BandWidth=parameter(1);    %����
Tpulse=parameter(2);       %ʱ��
Nazimuth=parameter(3);     %��λ�������
fs = parameter(4);         %��Ƶ
cj=sqrt(-1);
Kr=BandWidth/Tpulse;
T=Tpulse;
t1=linspace(-T/2,T/2,T*fs);
%%
% Filter=exp(-cj*pi*Kr*t1.^2);                  % ʱ�����˲���
N = length(t1);
Rwind =hamming(N);
% Filter1=exp(-cj*pi*Kr*t1.^2).*Rwind.';        	% ʱ�����˲���
Filter1=exp(-cj*pi*Kr*t1.^2);
TTT=Nrange-length(Filter1);
Filter=[Filter1 zeros(1,TTT)];
dataRR1 = fft(dataRR).*fft(Filter.'*ones(1,Nazimuth));     % Ƶ���˲�
Rwinp=([ones(1,4500) zeros(1,Nrange-9000) ones(1,4500)]);
dataRR2=dataRR1.*(Rwinp.'*ones(1,Nazimuth));
dataRR3 = ifft(dataRR2);

dataout=dataRR3;
end

