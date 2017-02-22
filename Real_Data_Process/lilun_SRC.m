function [ Srancom1 ] = lilun_SRC( varargin )
%  理论SRC
%  输入：parameter--由SRC_lilun_Parameter函数输出的参数列表
%        dataRR11--距离压缩后的回波数据 
%        SRC_lilun_Parameter=[a,Nrange,fc,SapTime];
% 输出：理论SRC
%% 参数读入
parameter=varargin{1};   %参数列表
Srancom=varargin{2};        %回波数据
Rrange = varargin{3};
fa = varargin{4};
fr = varargin{5};
a=parameter(1);  %带宽
RanNum=parameter(2);          %时宽
fc=parameter(3);
SapTime=parameter(4);
cj=sqrt(-1);
c=3e8;
%% SRC
load('hT.mat')
load('xT2.mat')
load('yT2.mat')
load('hR.mat')
load('xR2.mat')
load('yR2.mat')

Vr=[0 31 0];Vt=[0 31 0];
Location_T=[xT2(1) yT2(1)-(hR(1)-640+130)/tand(25) hT(1)-640+130];
Location_R=[xR2(1) yR2(1)-(hR(1)-640+130)/tand(25) hR(1)-640+130];

Rrange1=Rrange(a-999:a+1000);
x_T=Location_T(1);
y_T=Location_T(2);
z_T=Location_T(3);
x_R=Location_R(1);
y_R=Location_R(2);
z_R=Location_R(3);
R_bi=Rrange1;
options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off','TolFun',1e-9);
[x,F] = fsolve(@(x)sqrt((x-x_T).^2+y_T.^2+z_T.^2)+...
    sqrt((x-x_R).^2+y_R.^2+z_R.^2)-R_bi,zeros(1,RanNum),options);  %%场景中心点x=0，y=0；

TargetCenter=[x(RanNum/2+1) 0 0];
Rsq_T=sqrt((TargetCenter-Location_T)*(TargetCenter-Location_T).');
Rsq_R=sqrt((TargetCenter-Location_R)*(TargetCenter-Location_R).');
R_center=Rsq_T+Rsq_R;
Tr=SapTime*(-RanNum/2:RanNum/2-1);
theta_Rsq=asin(Vr*(TargetCenter-Location_R).'/sqrt(Vr*Vr.')/...
    sqrt((TargetCenter-Location_R)*(TargetCenter-Location_R).'));
theta_Tsq=asin(Vr*(TargetCenter-Location_T).'/sqrt(Vt*Vt.')/...
    sqrt((TargetCenter-Location_T)*(TargetCenter-Location_T).'));
thetaT=asin(Vr*(TargetCenter-Location_T).'/sqrt(Vt*Vt.')/...
    sqrt((TargetCenter-Location_T)*(TargetCenter-Location_T).'));
thetaR=asin(Vr*(TargetCenter-Location_R).'/sqrt(Vr*Vr.')/...
    sqrt((TargetCenter-Location_R)*(TargetCenter-Location_R).'));
Referent_Rsq_T=norm(TargetCenter-Location_T);
Referent_Rsq_R=norm(TargetCenter-Location_R);



k10=-Vr(2)*sin(thetaR)-Vt(2)*sin(thetaT);
k20=1/2*(Vr(2).^2*cos(thetaR).^2/Referent_Rsq_R+Vt(2).^2*cos(thetaT).^2/Referent_Rsq_T);
k30=1/6*(3*Vr(2).^3*cos(thetaR).^2*sin(thetaR)/Referent_Rsq_R.^2+...
    3*Vt(2).^3*cos(thetaT).^2*sin(thetaT)/Referent_Rsq_T.^2);
k40=1/24*(3*Vr(2).^4*cos(thetaR).^2*(5*sin(thetaR).^2-1)/Referent_Rsq_R.^3+...
    3*Vt(2).^4*cos(thetaT).^2*(5*sin(thetaT).^2-1)/Referent_Rsq_T.^3);
[Fr,Fa] = meshgrid(fr,fa);
phi_src = 2*pi*(c/4/k20/fc*((Fr/fc).^2-(Fr/fc).^3).*Fa.^2+...
    k30/8/k20^3*(3*k10*(c/fc)*((Fr/fc).^2-(Fr/fc).^3).*Fa.^2+...
    (c/fc)^2*(3*(Fr/fc).^2-4*(Fr/fc).^3).*Fa.^3)+(9*k30^2-...
    4*k20*k40)/64/k20^5*(6*k10^2*(c/fc)*((Fr/fc).^2-(Fr/fc).^3).*Fa.^2+...
    4*k10*(c/fc)^2*(3*(Fr/fc).^2-4*(Fr/fc).^3).*Fa.^3+...
    (c/fc)^3*(6*(Fr/fc).^2-10*(Fr/fc).^3).*Fa.^4));

Srancom1=Srancom.*exp(-cj*phi_src);
end

