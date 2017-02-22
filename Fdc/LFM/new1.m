clc;clear;close all;
%% --------------------------------帧起始32w----------------------------------------
%% 
%%%%%%%%%%%%%%%%%%%%    数据格式 RanNum * AziNum     %%%%%%%%%%%%%%%%%%%%%%%%%%
[ADpara.AD_filename,ADpara.AD_pathname]=uigetfile('E:\','\*.dat','请选择需要分析的AD数据');
fid=fopen(strcat(ADpara.AD_pathname,ADpara.AD_filename),'r');
% [ADParameter,DataOrg]=ReadAdFrameData0(fid);

n = 554;
AziNum = n*16;
RanNum = 8192;
echocell = zeros(8192,16);
echo = zeros(8192,AziNum);
for i=1:n
    [ADParameter,DataOrg]=ReadAdFrameData0(fid);
    thera_r11(i)=ADParameter.AnAz;
    thera_r22(i)=ADParameter.SAz;
    long(i)=ADParameter.pl_long;                                %经度
    lat(i)=ADParameter.pl_lat;                                  %纬度
    h(i)=ADParameter.pl_h;                                      %高度
    a=DataOrg{1,1};
    echocell=reshape(a,8192,16);
    echo(:,16*i-15:16*i) = echocell;
end
clear n i;
% figure();plot(h)
% max(h)-min(h)

F_echo=ftx(fty(echo));
figure('Name','回波');imagesc(abs(echo));colormap(gray)
figure('Name','回波频谱');imagesc(abs(F_echo));
figure;plot(abs((F_echo(4096,:))));colormap(gray)
%%
%%%%%%%%%%%%%%%%%%%%%%%    去多普勒质心    %%%%%%%%%%%%%%%%%
PRI=ADParameter.PRI*1e-6;
PRF=1/PRI;
Ta = PRI*(-AziNum/2:AziNum/2-1);  
% load fdcc.mat
fdccc=6836.5;
dataRR1=echo.*exp(-1i*2*pi*fdccc*ones(RanNum,1)*Ta);
figure('Name','回波');imagesc(abs(fty(ftx(dataRR1))));
F_echo=fty(ftx(dataRR1));
%% 
%%%%%%%%%%%%%%%%%%%%%%%   脉冲压缩      %%%%%%%%%%%%%%%%%%%%
Fs = 1/(ADParameter.Rw*1E-6);
fs=Fs;
Tr = ADParameter.Pw*1E-6;
Br = ADParameter.Bw*1E6;
Kr = Br/Tr;
dt = 1/Fs;


% RanTime = (-RanNum/2:RanNum/2-1)*dt;
% temp = abs(RanTime) < Tr/2;
% Size_All = size(temp,2);
% Size_1 = size(find(temp==1),2);
% Size_0 = Size_All - Size_1;
% M = (Size_1 - 1) / 2;
% k = -M:M;
% % win = 0.54-0.46*cos(2*pi*k/M); 
% win = ones(1,2*M+1);
% kall=[zeros(1,(Size_0-1)/2),win,zeros(1,(Size_0+1)/2)];
% windows = temp.*kall;
% F_MatchFilter = (exp(1i*pi*Kr*RanTime.^2).*windows).'*ones(1,AziNum);
% % figure('Name','匹配滤波');plot(real(s_ref(:,1)));
% s_out1 = (iftx((ftx(echo)).*(conj(ftx(F_MatchFilter)))));
% figure('Name','脉冲压缩');imagesc(abs(s_out1));colormap(gray);
% set(gca,'ydir','normal');


RanTime = (-RanNum/2:RanNum/2-1)*dt;
fr=Fs/RanNum*(-RanNum/2:RanNum/2-1);
F_MatchFilter = (exp(1i*pi/Kr*fr.^2)).'*ones(1,AziNum);
s_out=iftx(ifty(F_echo.*F_MatchFilter));
figure('Name','脉冲压缩');imagesc(abs((s_out)))
% s_out1=s_out(4400:5200,2500:6500);
%%
%%%%%%%%%%%%%%%%%%%%% radon变换估计fdc   %%%%%%%%%%%%%%%%%%%%%%%%
% % % % % 

% theta=0:179;
% R=radon(abs((s_out1)),theta);  %%%% R = radon(I, theta) 返回亮度图像在角度theta下的Radon变换R。
% % imagesc((R));          %%%%%估计出目标轨迹的斜率
% b=max(max((R)));
% [u,cc]=find((R)==b);
% cc=cc-1;
% theta=linspace(cc-1,cc+1,201);
% R=radon(abs((s_out1)),theta);
% b=max(max((R)));
% [u,c0]=find((R)==b);
% cc=cc-1+(c0-1)*0.01;     %%%%% c为目标真实方向的旋转角度
% if cc<=90
%     k=tan((90-cc)/180*pi);     %%%% k为估计出的目标斜率
% else
%     k=tan((270-cc)/180*pi);
% end
% c=3e8;
% lamda=c/9.6e9;
% V=ADParameter.pl_Ag;
% fdcc=-V*k*dt*c/PRI/V/lamda;     %%% 
% save('fdcc.mat','fdcc');
%%
%%%%%%%%%%%%%%%%%%%%%%%          RCMC                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRF=1/(ADParameter.PRI*1e-6);
% Srancom=fty(ftx(s_out));
% c=3e8;
% fs=Fs;
% SapTime=dt;
% delta_R = SapTime*c;
% K = tan(10*pi/180)*delta_R/(PRF/AziNum);    
% K1=tan(0.0*pi/180)*delta_R/(PRF/AziNum);
% t0 =PRF/AziNum*K*[-AziNum/2:AziNum/2-1]/c+K1/c*(PRF/AziNum*[-AziNum/2:AziNum/2-1]).^2;
% Srancom_jiaozheng=iftx((Srancom).*(exp(1i*pi*2*t0.'*(fs/RanNum*[-RanNum/2:RanNum/2-1]))).');
% Sran=ifty(Srancom_jiaozheng);
% figure;imagesc((abs((Sran))))


% KK=-0.1091465;
KK=-0.1146;
fr=fs/RanNum*[-RanNum/2:RanNum/2-1].';
coe=([1:AziNum])*(KK)*(1/fs);
data1=zeros(RanNum,AziNum);
for jjj=1:AziNum
    data1(:,jjj)=iftx(ftx(s_out(:,jjj)).*exp(1i*2*pi*coe(jjj)*fr));
end
figure('Name','RCMC');
imagesc(abs(data1)) 

figure;plot(imag(data1(2090,:)))
figure;plot(real(data1(2090,:)))
% 
% figure;plot(unwrap(angle(data1(2090,:))))
%%
%%%%%%%%%%%%%%%%%%           估计调频斜率               %%%%%%%%%%%%%%%%
x1=data1(2090,3401:4800);
Na=1400;
Tsar=Na/PRF;
[C,D]=tfrate(x1,1/PRF,Tsar); %对消后动目标所在距离单元的时间—调频斜率分布 
G=abs(C);
figure
colormap(gray(256))
xg=max(max(G));ng=min(min(G));cg=255/(xg-ng);
imagesc(cg*(G-ng))                   %显示时间—调频斜率分布图     
figure
plot(10*log10(D/max(D)))                 %画出时间—调频斜率沿时间轴积分的结果
xlabel('调频率')
ylabel('归一化幅度')
zoom xon;
v=find(abs(D)>=max(abs(D)));
fdr=(2/PRF/Tsar/Na)*(v-Na/2-1)*PRF^2; 
%%
%%%%%%%%%%%%%%%%%%%%%                方位压缩                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:10
%     i
%     fdr=-40-i*1;
% Ta=(-AziNum/2:AziNum/2-1)*(1/PRF);
% Refa1=(((exp(1i*pi*(fdr)*Ta.^2)))).';
% A=fty(data1);
% B=(conj(ftx(Refa1))*ones(1,RanNum));
% Imageraw=ifty(A.*(B.'));
% e(i)=my_entropy(Imageraw);
% % figure;imagesc(abs((Imageraw)))
% end

Ta=(-AziNum/2:AziNum/2-1)*(1/PRF);
% for i=1:20
%     i
% fdr=-14+(i)*0.2;

Refa1=(((exp(1i*pi*(-42)*Ta.^2-1i*pi*0*Ta)))).';
A=fty(data1);
B=(conj(ftx(Refa1))*ones(1,RanNum));
Imageraw=ifty(A.*(B.'));
% e(i)=my_entropy(Imageraw);
% end
% figure;plot(e);
figure('Name','方位压缩');imagesc(abs((Imageraw)));colormap(gray)



Refa1=data1(2090,:).';
A=fty(data1);
B=(conj(ftx(Refa1))*ones(1,RanNum));
Imageraw=ifty(A.*(B.'));
% e(i)=my_entropy(Imageraw);
% end
% figure;plot(e);
figure('Name','方位压缩');imagesc(abs((Imageraw)));colormap(gray)
%%
%%%%%%%%%%%%%%%%%     ps存图 .raw       %%%%%%%%%%%%%%
pixel_max=max(max(abs(Imageraw)));
pixel_min=min(min(abs(Imageraw)));
image=(abs(Imageraw)-pixel_min)*255/(pixel_max-pixel_min);
uplimite=255;
image=image.*(image<uplimite)+uplimite*(image>=uplimite);
pixel_max=max(max(abs(image)));
pixel_min=min(min(abs(image)));
image=(abs(image)-pixel_min)*255/(pixel_max-pixel_min);

fid=fopen('p1.raw','wb');
fwrite(fid,uint8(image.'),'uint8');
fclose('all');

%%
ffs=9.6e9;
c=3e8;
lambda=c/ffs;
R0=ADParameter.R0*c*1e-6;
V=ADParameter.pl_Ag;
theta_r1=ADParameter.AnAz;
theta_r2=ADParameter.SAz;
fdr=-2*V*V*(cos(pi/2-theta_r1))^2/lambda/R0;
fdc=2*V*sin(pi/2-theta_r1)/lambda;


