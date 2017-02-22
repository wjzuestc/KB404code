function [dataout,echoData] = ReaddataFMF( varargin)
%  读取回波函数
%  输入：parameter--由Readparameter函数输出的参数列表
%           dijiao--待读取的第dijisao扫
%          Nrange1--距离向采样点数（包括帧头）
%           select--用于确定第dijisao扫是否和场景成镜像关系
% parameter=[fuyangjiao;Bandwidth;T;PRF;Vscan;Amin;Amax;H;Va;AngleVH;AH;AV;AR;Longitude;Latitude;dushift;Nazimuth];
parameter=varargin{1};
read_shift=varargin{2};
dijisao=varargin{3};
Nrange1=varargin{4};
select=varargin{5};
dataPath = varargin{6};
switch num2str(nargin)
    case '7'
        sample_beishu=varargin{7};
    case '8'
        Nmin=varargin{7};
        Nmax=varargin{8};
    case '9'
        Nmin=varargin{7};
        Nmax=varargin{8};
        sample_beishu=varargin{9};
    otherwise
        error('输入参数不对')
end
Bandwidth=parameter(2);
T=parameter(3);
Nazimuth=parameter(5);
Nblock=parameter(8);
samplerate_full = parameter(7);
fc = parameter(6);           % 载频
c=3e8;
cj=sqrt(-1);
%%%%%%%%%%%%%%% 判断是否满足实信号带通采样定理
fH=fc+Bandwidth*1.2/2;
fL=fc-Bandwidth*1.2/2;
fs=samplerate_full/sample_beishu; % 采样率
Ms=linspace(0,floor(fH/Bandwidth/1.2),floor(fH/Bandwidth/1.2)+1);
fsH=2*fL./Ms;
fsL=2*fH./(Ms+1);
count=0;
for i=1:length(Ms)
    if fsL(i)<fs&&fs<fsH(i)
        count=count+1;
    end
end
if count==0
    error('不满足实信号带通采样定理');
end
%%%%%%%%%%%%%%%
i=dijisao;
Binhead=floor((read_shift+(i-1)*Nazimuth)/Nblock);
if rem((read_shift+i*Nazimuth),Nblock)==0
    Bintail=(read_shift+i*Nazimuth)/Nblock-1;
else
    Bintail=floor((read_shift+i*Nazimuth)/Nblock);
end
if Bintail-Binhead==0
    name=[dataPath,'\',num2str(Binhead),'.bin'];
    handle = fopen(name);
    x=read_shift+Nazimuth*(i-1)-Binhead*Nblock;
    fseek(handle,x*Nrange1,'bof');
    data1=fread(handle,[Nrange1,Nazimuth],'int8');
    fclose(handle);
    %######################
    %读帧头
    %######################
    frame=zeros(6,Nazimuth);
    handle = fopen(name);
    for xx=1:Nazimuth
        fseek(handle,x*Nrange1+(xx-1)*Nrange1,'bof');
        frame(:,xx)=fread(handle,6,'uint8');
    end
    fclose(handle);
else if Bintail-Binhead==1
        name=[dataPath,'\',num2str(Binhead),'.bin'];
        handle = fopen(name);
        x=read_shift+Nazimuth*(i-1)-Binhead*Nblock;
        fseek(handle,x*Nrange1,'bof');
        data2=fread(handle,[Nrange1,Nblock-x],'int8');
        fclose(handle);
        %######################
        %读帧头
        %######################
        handle = fopen(name);
        frame2=zeros(6,Nblock-x);
        for xx=1:Nblock-x
            fseek(handle,x*Nrange1+(xx-1)*(Nrange1),'bof');
            frame2(:,xx)=fread(handle,6,'uint8');
        end
        fclose(handle);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        y=read_shift+Nazimuth*(i)-Bintail*Nblock;
        name=[dataPath,'\',num2str(Bintail),'.bin'];
        handle = fopen(name);
        fseek(handle,0,'bof');
        data3=fread(handle,[Nrange1,y],'int8');
        fclose(handle);
        data1=[(data2),(data3)];
        handle = fopen(name);
        frame3=zeros(6,y);
        for xx=1:y
            fseek(handle,(xx-1)*Nrange1,'bof');
            frame3(:,xx)=fread(handle,6,'uint8');
        end
        fclose(handle);
        frame=[frame2,frame3];
    else if Bintail-Binhead==2
            name=[dataPath,'\',num2str(Binhead),'.bin'];
            handle = fopen(name);
            x=read_shift+Nazimuth*(i-1)-Binhead*Nblock
            fseek(handle,x*Nrange1,'bof');
            data2=fread(handle,[Nrange1,Nblock-x],'int8');
            fclose(handle);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            name=[dataPath,'\',num2str(Binhead+1),'.bin'];
            handle = fopen(name);
            fseek(handle,0,'bof');
            data3=fread(handle,[Nrange1,Nblock],'int8');
            fclose(handle);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y=read_shift+Nazimuth*(i)-Bintail*Nblock
            fseek(handle,0,'bof');
            name=[dataPath,'\',num2str(Bintail),'.bin'];
            handle = fopen(name);
            data4=fread(handle,[Nrange1,y],'int8');
            data1=[data2,data3,data4];
            fclose(handle);
            frame=zeros(6,Nazimuth);
            for xx=1:Nazimuth
                fseek(handle,x*Nrange1+(xx-1)*Nrange1,'bof');
                frame(:,xx)=fread(handle,6,'uint8');
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%555
% data1 =data2;
frame1=zeros(1,Nazimuth);
for xx=1:Nazimuth
    frame1(1,xx)=frame(6,xx)*2^24+frame(5,xx)*2^16+frame(4,xx)*2^8+frame(3,xx);
end
check=frame1(Nazimuth)-frame1(1)+1-Nazimuth
data=data1(101:Nrange1,:);
[Mdata,Ndata]=size(data);
Nrange=Mdata;
Nazimuth=Ndata;
data1=downsample(data,sample_beishu);
data=data1;
[M,N]=size(data);
Nrange=M;
% hangzi = zeros(32,1);
% fileTest = fopen('test.bin','wb');
% for kk = 1:N
%     fwrite(fileTest,frame2(:,kk),'uint8');
%     fwrite(fileTest,hangzi,'uint8');
%     fwrite(fileTest,data(:,kk),'int8');
% end
% fclose(fileTest);
%######################
%希尔伯特变换
%######################
Kr=Bandwidth/T;
data = (fft(data));
basefc=rem(fc,fs);
if basefc<fs/2
    Filter = [zeros(1,ceil(Nrange/2)+100),ones(1,floor(Nrange/2)-100)];
    fc_shift=basefc;
else
    Filter = [ones(1,ceil(Nrange/2)+100),zeros(1,floor(Nrange/2)-100)];
    %         Filter = [zeros(1,ceil(Nrange/2)+100),ones(1,floor(Nrange/2)-100)];
    fc_shift=-fs+basefc;
end
hh = fftshift(Filter);
data = data.*(fftshift(Filter)'*ones(1,Nazimuth));
data = ifft((data));
echoData=data;
%######################
%去载频
%######################
t=(0:Nrange-1)/fs; % 时间轴，用于构造滤波器
data = data.*(exp(-cj*2*pi*fc_shift*t).'*ones(1,Nazimuth));
% end
if nargin==8
    dataout=data(Nmin:Nmax,:);
else
    dataout=data;
end
end

