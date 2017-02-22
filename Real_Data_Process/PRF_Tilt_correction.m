function  S1  = PRF_Tilt_correction( varargin )
%  PRF倾斜校正
%  输入：parameter--由freRCMC_Parameter函数输出的参数列表
%        dataRR11--距离压缩后的回波数据 
%        freRCMC_Parameter = [AziNum,fr,RanNum,fs,SapTime];
% 输出：PRF校正后的数据
%% 参数读入
parameter=varargin{1};   %参数列表
S=varargin{2};    %回波数据
fr=varargin{3};
c=3e8;
cj=sqrt(-1);
AziNum=parameter(1);
RanNum=parameter(2);
%% PRF倾斜校正
Fr=ones(AziNum,1)*fr;
Rerror = linspace(0,123.84,AziNum).'*ones(1,RanNum);        %%实际情况和理想情况到场景中心的距离历史和之差     
S1 =ifty(fty(S).*exp(cj*2*pi/c*Rerror.*Fr));        %%距离向重采样（解决回波存储位置错位跨越距离单元）
% S2=S1(:,1001:4000);
end

