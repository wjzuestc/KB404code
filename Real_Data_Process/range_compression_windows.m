function [ dataout ] = range_compression_windows( varargin )
%  对回波进行距离向脉冲压缩--加窗
%  输入：parameter--由Readparameter函数输出的参数列表
%        data--回波数据 
%        Nrange=varargin{3};
% parameter=[fuyangjiao;Bandwidth;T;PRF;Vscan;Amin;Amax;H;Va;AngleVH;AH;AV;AR;Longitude;Latitude;dushift;Nazimuth];
% 输出：距离向脉冲压缩（加窗）后的回波数据
%% 参数收入
parameter=varargin{1};     %参数列表
dataRR=varargin{2};        %回波数据
Nrange=varargin{3};        %距离向采样点数
BandWidth=parameter(1);    %带宽
Tpulse=parameter(2);       %时宽
Nazimuth=parameter(3);     %方位向采样点
fs = parameter(4);         %载频
cj=sqrt(-1);
Kr=BandWidth/Tpulse;
T=Tpulse;
t1=linspace(-T/2,T/2,T*fs);
%%
% Filter=exp(-cj*pi*Kr*t1.^2);                  % 时域构造滤波器
N = length(t1);
Rwind =hamming(N);
% Filter1=exp(-cj*pi*Kr*t1.^2).*Rwind.';        	% 时域构造滤波器
Filter1=exp(-cj*pi*Kr*t1.^2);
TTT=Nrange-length(Filter1);
Filter=[Filter1 zeros(1,TTT)];
dataRR1 = fft(dataRR).*fft(Filter.'*ones(1,Nazimuth));     % 频域滤波
Rwinp=([ones(1,4500) zeros(1,Nrange-9000) ones(1,4500)]);
dataRR2=dataRR1.*(Rwinp.'*ones(1,Nazimuth));
dataRR3 = ifft(dataRR2);

dataout=dataRR3;
end

