function [ Imageraw ] = azimuth_compression( varargin )
%  ��λ������ѹ��
%  ���룺parameter--��azimuth_compression_Parameter��������Ĳ����б�
%        dataRR11--����ѹ����Ļز����� 
%        azimuth_compression_Parameter = [fdr,Ta,Ts,Nrange];
% �������λ������ѹ���������
%%
fdr=varargin{1};   %�����б�
Ta=varargin{2}; 
Ts=varargin{3}; 
Nrange=varargin{4}; 
S=varargin{5};    %�ز�����
cj=sqrt(-1);
%%
%��λѹ��
A=((exp(cj*pi*(fdr+45)*Ta.^2-cj*pi*(0)*Ta)));
B=( abs(Ta)<Ts/2 );
Refa1=(A.*B).';
A=ftx(S);
B=(conj(ftx(Refa1))*ones(1,Nrange));
Imageraw=A.*B;
figure;colormap(gray);imagesc(abs((Imageraw.')))
end

