function [dataout] = range_compression( varargin )
%  对回波进行距离向脉冲压缩（不加窗）
%  输入：parameter--由Readparameter函数输出的参数列表
%        data--回波数据 
%        range_compression_Parameter=[BandWidth;Tpulse;AziNum;fs];
% 输出：距离向脉冲压缩后的回波数据
%% 参数收入
parameter=varargin{1};   %参数列表
data=varargin{2};        %回波数据
Bandwidth=parameter(1);  %带宽
T=parameter(2);          %时宽
Nazimuth=parameter(3);
samplerate_full = parameter(4);

sample_beishu = 1;   %？？？
cj=sqrt(-1);
fs=samplerate_full/sample_beishu; % 采样率
Kr=Bandwidth/T;                   %调频斜率
[Mdata,Ndata]=size(data); 
Nrange=Mdata;
%% 距离向脉冲压缩不加窗
t1=linspace(-T/2,T/2,T*fs);
Filter=exp(-cj*pi*Kr*t1.^2);                  % 时域构造滤波器
TTT=Nrange-length(Filter);
Filter=[Filter zeros(1,TTT)];
data = fft(data).*fft(Filter.'*ones(1,Nazimuth));     % 频域滤波
data = ifft(data);

dataout=data;
end

