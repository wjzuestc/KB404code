function [ Srancom_jiaozheng ] = time_domain_RCMC( varargin )
%  时域进行RCMC
%  输入：parameter--由timeRCMC_Parameter函数输出的参数列表
%        dataRR11--距离压缩后的回波数据 
%        timeRCMC_Parameter = [fdc,lambda,V,Ta,fr,fdr];
% 输出：时域RCMC后的数据
%% 参数收入
parameter=varargin{1};   %参数列表
dataRR11=varargin{2};        %回波数据
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
K0=fdc/1.5*lambda/V;%%%%距离徙动校正需要调试的参数%%%%%%%%
t_ij=-V*(K0)*Ta/c;
comf=exp(cj*2*pi*t_ij.'*fr);             
Sy1=Sy.*comf;
S1=ifty(Sy1);


K1=fdr/8*lambda;%%%%距离徙动校正需要调试的参数%%%%%%%%
t_ij1=-(K1)*Ta.^2/c;
comf=exp(cj*2*pi*t_ij1.'*fr);             
Sy2=Sy1.*comf;
Srancom_jiaozheng=ifty(Sy2);
end

