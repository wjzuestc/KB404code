function [ Srancom1 ] = SRC_estimation( varargin )
%  估计SRC
%  输入：parameter--由SRC_estimation_Parameter函数输出的参数列表
%        dataRR11--距离压缩后的回波数据 
%        SRC_estimation_Parameter=[BandWidth,Tpulse,AziNum];
% 输出：SRC处理后的数据
%% 参数读入
parameter=varargin{1};   %参数列表
Srancom=varargin{2};        %回波数据
fr=varargin{3};
BandWidth=parameter(1);  %带宽
Tpulse=parameter(2);          %时宽
AziNum=parameter(3);
Kr=BandWidth/Tpulse;
%%
Filter=exp(j*pi*fr.^2/Kr/6);                  % 频域构造滤波器
Srancom1=Srancom.*(ones(AziNum,1)*Filter);
end

