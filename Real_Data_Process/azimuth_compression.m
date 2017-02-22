function [ Imageraw ] = azimuth_compression( varargin )
%  方位向脉冲压缩
%  输入：parameter--由azimuth_compression_Parameter函数输出的参数列表
%        dataRR11--距离压缩后的回波数据 
%        azimuth_compression_Parameter = [fdr,Ta,Ts,Nrange];
% 输出：方位向脉冲压缩后的数据
%%
fdr=varargin{1};   %参数列表
Ta=varargin{2}; 
Ts=varargin{3}; 
Nrange=varargin{4}; 
S=varargin{5};    %回波数据
cj=sqrt(-1);
%%
%方位压缩
A=((exp(cj*pi*(fdr+45)*Ta.^2-cj*pi*(0)*Ta)));
B=( abs(Ta)<Ts/2 );
Refa1=(A.*B).';
A=ftx(S);
B=(conj(ftx(Refa1))*ones(1,Nrange));
Imageraw=A.*B;
figure;colormap(gray);imagesc(abs((Imageraw.')))
end

