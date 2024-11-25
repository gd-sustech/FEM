%%%% Copyright:
%%%% You can copy and use this code for any purpose, but you must declare 
%%%% that the original version of this code was written by zxf
%%%% RB fluid computed with mrtLBM  
% from ... Soret and Dufour
%add square solids
%------------------------------------------------------------------------------------------------
%% Programming Initialization
close all; clear; clc;
fprintf('Programming Rayleigh Benard with tri-Mrt model Running ... \n');
%load porous.04.mat
%load porous.036H1_full.mat

Nx=1001;Ny=1001; Lx=1;Ly=1;
F=zeros(Nx,Ny,2);%欧拉点受力情况

delta_x=Lx/(Nx-1);delta_y=delta_x; 
x = (0:(Nx-1))*delta_x;
y = (0:(Ny-1))*delta_y;
[X,Y] = meshgrid(x,y);


G=zeros(Nx,Ny);
Tf=zeros(Nx,Ny);   
%lamta=0;nerror=1;%ibm松弛系数及迭代次数
lamta=2.5;nerror=10;Ferror=zeros(Nx,Ny,2,nerror);
Terror=zeros(Nx,Ny,nerror);
% Model's parameters
Dimension = 2;                    % Problem's Dimension
Q = 9;                            % Qian 9 D2Q9 Velocity model 
%U = 0.05;                          % velocity of inlet
U = 0;
rho_0 = 1.;                       % reference density

% Discrete velocity :
ev = [0  1  0 -1  0  1 -1 -1  1
     0  0  1  0 -1  1  1 -1 -1 ]; 
evS=[0 1 0 -1 0;0 0 1 0 -1];
% weight in equilibrium distribution function :
omega= [4./9. 1./9. 1./9. 1./9. 1./9. ...
         1./36. 1./36. 1./36. 1./36.];  
%omegaS=[-1/5 3/10 3/10 3/10 3/10];
    %omegaS=[3/5 1/10 1/10 1/10 1/10];   
    omegaS=[1/3 1/6 1/6 1/6 1/6];   
    %A multiple-relaxation-time lattice Boltzmann model for...
    ...convection heat transfer in porous media
c = 1 ;                                     % speed of lattice sound    (must equal to 1)
delta_t = delta_x / c ;                       % time step interval
%mu=0.01;
Pr=5.43;

%Ra=10^8*9.81;N=1.5; %验证case
%Ra=2*10^8*9.81*100;N=1.5; %验证case

%mu=2e-3;Ra=2e6;
mu=5e-5;
Ra=1e9;

c_square = c^2. ;
a1 = 3. / c_square ; 
a2 = 9./2. / c_square^2. ;
a3 =  3./2. / c_square ;
tau_f = 3. * mu / (c^2. * delta_t) + 1./2. ;
tau_h=(3. * mu/Pr / (c^2. * delta_t) + 1./2.);
gbeta=Ra*mu*(mu/Pr)/Ly^3; 

% Macroscopic 
rho = zeros(Nx,Ny);                       % macro density
u = zeros(Nx,Ny,2);                       % macro velocity
u_temp = u;                                   % macro temporary velocity

feq = zeros(Q,1);                             % equilibrium distribution function
tmax=200000;
% Outputs control
%OutputsInterval = 1;
OutputsInterval = 500;
files_name = 'salt finger';
% Plot control
% PlotInterval = 1;
%PlotInterval = 50;
PlotInterval = 100;
plot_name1 = 'Velocity Field'; 
plot_name2 = 'Pressure Field';
feq1=zeros(Nx-1,Ny-1);
edotu1=zeros(Nx-1,Ny-1);
usq=zeros(Nx+1,Ny+1);
ii=2:Nx-1;jj=2:Ny-1;

%    shizongx=zeros(200,tmax+1);shizongy=zeros(200,tmax+1);
%    shizongx(:,1)=0.01:0.01:2;shizongy(:,1)=repmat(((Ny+1)/2)*delta_y,200,1);

%四个方向，Flag变量，Plot变量
%设计多孔介质固体骨架：
ssw=0.036;ssh=0.036; SSW=ssw/delta_x+1;SSH=ssh/delta_x+1;
scx=0.025:0.05:1;
scy=0.025:0.05:1;
[Scx,Scy]=ndgrid(scx,scy);
sxb=Scx-ssw/2;  SXB=sxb/delta_x+1; SXB=SXB(:);
sxe=Scx+ssw/2;  SXE=sxe/delta_x+1; SXE=SXE(:);
syb=Scy-ssh/2;  SYB=syb/delta_x+1; SYB=SYB(:);
sye=Scy+ssh/2;  SYE=sye/delta_x+1; SYE=SYE(:);

nS=length(Scx(:));
locb=zeros(1,nS*SSW);
loct=zeros(1,nS*SSW);
locl=zeros(1,nS*SSH);
locr=zeros(1,nS*SSH);

Flag_squ=zeros(Nx,Ny);

for flag=1:length(Scx(:))
    locb((1:SSW)+(flag-1)*SSW)=(SYB(flag)-1)*Nx+(round(SXB(flag)):round(SXE(flag))); 
    loct((1:SSW)+(flag-1)*SSW)=(SYE(flag)-1)*Nx+(round(SXB(flag)):round(SXE(flag)));
    locl((1:SSH)+(flag-1)*SSH)=(round(SYB(flag)):round(SYE(flag)))*Nx+SXB(flag)-Nx;
    locr((1:SSH)+(flag-1)*SSH)=(round(SYB(flag)):round(SYE(flag)))*Nx+SXE(flag)-Nx;
    Flag_squ(X'>=sxb(flag)&X'<=sxe(flag)&Y'>=syb(flag)&Y'<=sye(flag))=1;
end
locb=round(locb); loct=round(loct); locl=round(locl); locr=round(locr);

%% Initial the flow field
rho=ones(Nx,Ny)*rho_0;
u(:,:,:) = 0;
T=zeros(Nx,Ny);
%T=0.5*ones(Nx,Ny);

%{
%S(ii,((Ny-1)/4+10))=0.5+0.5*sin((1:800)/160*2*pi)+0.04*(rand(1,800)-0.5);
suiji=rand(1,800);
S(ii,((Ny-1)/2+5))=suiji;
S(1,((Ny-1)/2+5))=S(Nx-1,((Ny-1)/2+5));S(Nx,((Ny-1)/2+5))=S(2,((Ny-1)/2+5));
S(:,((Ny-1)/2+5)+1)=S(:,((Ny-1)/2+5)+1);S(:,((Ny-1)/2+5)-1)=S(:,((Ny-1)/2+5)-1);
%}
f = zeros(Nx,Ny,Q);
h=zeros(Nx,Ny,5);

for flag=1:9
    f(:,:,flag)=omega(flag)*rho;   
end


MT=[1 1 1 1 1;0 1 0 -1 0;0 0 1 0 -1;0 1 1 1 1;0 1 -1 1 -1];
%MT=[1 1 1 1 1;0 1 -1 0 0;0 0 0 1 -1;4 -1 -1 -1 -1;0 1 1 -1 -1];
%s1T=[0 1/(3*mu/Pr+0.5) 1/(3*mu/Pr+0.5) 1/(3^0.5/3+0.5) 1/(3^0.5/3+0.5)];
%s1T=[1 1/tau_h 1/tau_h 1.2 1.2];
s1T=[1 1/tau_h 1/tau_h 1.9 1.9];
MinvT=inv(MT);
Minv_sT=MinvT*diag(s1T);
%{
for flag=1:5
    s(:,:,flag)=omegaS(flag)*S;
    h(:,:,flag)=omegaS(flag)*T;
end
%}    
%

f1T=permute(h,[3,1,2]);
mT=MT*f1T(1:5,:);
     mT_eq(1,:,:)=T;
     mT_eq(2,:,:)=T.*u(:,:,1);
     mT_eq(3,:,:)=T.*u(:,:,2);
     %mT_eq(4,:,:)=-2*T/3;
     mT_eq(4,:,:)=2*T/3;   %zxf
     mT_eq(5,:,:)=0;
mT_temp0=MinvT*mT_eq(1:5,:);
mT_temp=reshape(mT_temp0,5,Nx,Ny);f1T=mT_temp;
h=permute(f1T,[2,3,1]);
%}


f_temp = f;h_temp=h;
M=[1 1 1 1 1 1 1 1 1;-4 -1 -1 -1 -1 2 2 2 2;4 -2 -2 -2 -2 1 1 1 1; ...
       0 1 0 -1 0 1 -1 -1 1;0 -2 0 2 0 1 -1 -1 1;0 0 1 0 -1 1 1 -1 -1; ...
       0 0 -2 0 2 1 1 -1 -1;0 1 -1 1 -1 0 0 0 0;0 0 0 0 0 1 -1 1 -1;];
s1=zeros(9,1);
         %s1(1)至s1(7)在（0，2）！
         %s1(1)=0; s1(4)=0;s1(6)=0;
         s1(1)=1; s1(4)=1;s1(6)=1;    
        % s1(2)=1.64; s1(3)=1.54; %no special restrictions
        s1(2)=1.9305;s1(3)=1.9305;   %no special restrictions 
        s1(2)=1.64;s1(3)=1.54;
         s1(8)=1/tau_f;  
         s1(9)=1/tau_f;
         %s1(5)=(8-s1(8))/8/(2-s1(8));  %for halfway bounce back
         %rule;porous media dsj=3/16
         %s1(5)=4*(2-s1(8))/(4+7*s1(8)); %for IBM dsj=9/8
         dsj=3/16;
         s1(5)=1/(dsj/(1/s1(8)-0.5)+0.5);
         %s1(5)=1.9;  %niu2006
         %s1(5)=s1(8);  %niu2006
         s1(7)=s1(5);
         %s1(:)=s1(8);%try single
        s2=repmat(s1,1,Nx*Ny);         
sss=diag(s1);%对角矩阵-松弛因子
Minv=inv(M);
Minv_s=Minv*sss;
m_eq=zeros(9,Nx,Ny);
Fi=zeros(9,Nx,Ny);
m=zeros(9,Nx*Ny);%zxf_les
m_temp0=zeros(9,Nx*Ny);
Feq=zeros(Nx,Ny,Q);


%tmax=1;

%{
writerObj=VideoWriter('RBinPorous_Ra4e8f51s07.avi');%1
open(writerObj);%1
set(gca,'nextplot','replacechildren');%1
set(gcf,'Renderer','zbuffer');%1
%}

Nu=zeros(tmax/10,1);
Re=zeros(tmax/10,1);
Txy=zeros(tmax/10,2);
L=zeros(tmax/10,1);
Area=sum(sum(Flag_squ == 0));

%for conjugate heat transfer


for t=1:tmax
t
%sum(sum(T))/Nx/Ny
%
%ceshi
%F=zeros(Nx,Ny,2);
%G=gbeta*((T-repmat(sum(T.*Flag_cir)./(1e-6+sum(Flag_cir)),Nx,1))-...
  %  N*(S-repmat(sum(S.*Flag_cir)./(1e-6+sum(Flag_cir)),Nx,1))).*rho;
%▲T
%G=gbeta*((T-repmat(mean(T),Nx,1))).*rho;
%G=gbeta*((T-0.5)).*rho;
G=gbeta*((T)).*rho;

F=zeros(Nx,Ny,2);
Fx=F(:,:,1);Fy=F(:,:,2)+G;



f1=permute(f,[3,1,2]);
m=M*f1(1:9,:);
     m_eq(1,:,:)=rho;
     m_eq(2,:,:)=rho.*(-2+3.*(u(:,:,1).^2+u(:,:,2).^2));
     %m_eq(3,:,:)=rho.*(1-2.*(u(:,:,1).^2+u(:,:,2).^2)); %by hu
     m_eq(3,:,:)=rho.*(1-3.*(u(:,:,1).^2+u(:,:,2).^2)); %by lu
     m_eq(4,:,:)=rho.*u(:,:,1);
     m_eq(5,:,:)=-rho.*u(:,:,1);
     m_eq(6,:,:)=rho.*u(:,:,2);
     m_eq(7,:,:)=-rho.*u(:,:,2);
     m_eq(8,:,:)=rho.*(u(:,:,1).^2-u(:,:,2).^2);
     m_eq(9,:,:)=rho.*u(:,:,1).*u(:,:,2);
     Fi(1,:,:)=0;
     Fi(2,:,:)=6*(Fx.*u(:,:,1)+Fy.*u(:,:,2));
     Fi(3,:,:)=-Fi(2,:,:);
     Fi(4,:,:)=Fx;
     Fi(5,:,:)=-Fx;
     Fi(6,:,:)=Fy;
     Fi(7,:,:)=-Fy;
     Fi(8,:,:)=2*(Fx.*u(:,:,1)-Fy.*u(:,:,2));
     Fi(9,:,:)=Fx.*u(:,:,2)+Fy.*u(:,:,1);
     m_temp0=s2.*(m-m_eq(1:9,:))-(1-s2/2)*delta_t.*Fi(1:9,:);
     m_temp0=Minv*m_temp0;
          m_temp=reshape(m_temp0,9,Nx,Ny);f1=f1-m_temp;
     %f1(:)=f1(:)-m_temp(:);   
f=permute(f1,[2,3,1]);
k=1;f(:,:,k) =circshift(f(:,:,k),[ev(1,k),ev(2,k),0]);
k=2;f(:,:,k) =circshift(f(:,:,k),[ev(1,k),ev(2,k),0]);
k=3;f(:,:,k) =circshift(f(:,:,k),[ev(1,k),ev(2,k),0]);
k=4;f(:,:,k) =circshift(f(:,:,k),[ev(1,k),ev(2,k),0]);
k=5;f(:,:,k) =circshift(f(:,:,k),[ev(1,k),ev(2,k),0]);
k=6;f(:,:,k) =circshift(f(:,:,k),[ev(1,k),ev(2,k),0]);
k=7;f(:,:,k) =circshift(f(:,:,k),[ev(1,k),ev(2,k),0]);
k=8;f(:,:,k) =circshift(f(:,:,k),[ev(1,k),ev(2,k),0]);
k=9;f(:,:,k) =circshift(f(:,:,k),[ev(1,k),ev(2,k),0]);
%3,求f
%collision+streaming

nn1=5;nn2=3;fff=f(:,:,nn1);fff=fff(:);fff2=f(:,:,nn2);fff2=fff2(:);fff(locb)=fff2(locb);f(:,:,nn1)=reshape(fff,Nx,Ny);
nn1=9;nn2=7;fff=f(:,:,nn1);fff=fff(:);fff2=f(:,:,nn2);fff2=fff2(:);fff(locb)=fff2(locb);f(:,:,nn1)=reshape(fff,Nx,Ny);
nn1=8;nn2=6;fff=f(:,:,nn1);fff=fff(:);fff2=f(:,:,nn2);fff2=fff2(:);fff(locb)=fff2(locb);f(:,:,nn1)=reshape(fff,Nx,Ny);
nn1=3;nn2=5;fff=f(:,:,nn1);fff=fff(:);fff2=f(:,:,nn2);fff2=fff2(:);fff(loct)=fff2(loct);f(:,:,nn1)=reshape(fff,Nx,Ny);
nn1=6;nn2=8;fff=f(:,:,nn1);fff=fff(:);fff2=f(:,:,nn2);fff2=fff2(:);fff(loct)=fff2(loct);f(:,:,nn1)=reshape(fff,Nx,Ny);
nn1=7;nn2=9;fff=f(:,:,nn1);fff=fff(:);fff2=f(:,:,nn2);fff2=fff2(:);fff(loct)=fff2(loct);f(:,:,nn1)=reshape(fff,Nx,Ny);
nn1=7;nn2=9;fff=f(:,:,nn1);fff=fff(:);fff2=f(:,:,nn2);fff2=fff2(:);fff(locl)=fff2(locl);f(:,:,nn1)=reshape(fff,Nx,Ny);
nn1=4;nn2=2;fff=f(:,:,nn1);fff=fff(:);fff2=f(:,:,nn2);fff2=fff2(:);fff(locl)=fff2(locl);f(:,:,nn1)=reshape(fff,Nx,Ny);
nn1=8;nn2=6;fff=f(:,:,nn1);fff=fff(:);fff2=f(:,:,nn2);fff2=fff2(:);fff(locl)=fff2(locl);f(:,:,nn1)=reshape(fff,Nx,Ny);
nn1=6;nn2=8;fff=f(:,:,nn1);fff=fff(:);fff2=f(:,:,nn2);fff2=fff2(:);fff(locr)=fff2(locr);f(:,:,nn1)=reshape(fff,Nx,Ny);
nn1=2;nn2=4;fff=f(:,:,nn1);fff=fff(:);fff2=f(:,:,nn2);fff2=fff2(:);fff(locr)=fff2(locr);f(:,:,nn1)=reshape(fff,Nx,Ny);
nn1=9;nn2=7;fff=f(:,:,nn1);fff=fff(:);fff2=f(:,:,nn2);fff2=fff2(:);fff(locr)=fff2(locr);f(:,:,nn1)=reshape(fff,Nx,Ny);

%4,f边界
%
rho(:,:)=sum(f,3);
u(:,:,1)=f(:,:,1)*ev(1,1)+f(:,:,2)*ev(1,2)+f(:,:,3)*ev(1,3)+f(:,:,4)*ev(1,4)+f(:,:,5)*ev(1,5)+...
        f(:,:,6)*ev(1,6)+f(:,:,7)*ev(1,7)+f(:,:,8)*ev(1,8)+f(:,:,9)*ev(1,9)+F(:,:,1)*delta_t/2;
u(:,:,2)=f(:,:,1)*ev(2,1)+f(:,:,2)*ev(2,2)+f(:,:,3)*ev(2,3)+f(:,:,4)*ev(2,4)+f(:,:,5)*ev(2,5)+...
        f(:,:,6)*ev(2,6)+f(:,:,7)*ev(2,7)+f(:,:,8)*ev(2,8)+f(:,:,9)*ev(2,9)+(F(:,:,2)+G)*delta_t/2;
u(:,:,1)=u(:,:,1)./rho(:,:);
u(:,:,2)=u(:,:,2)./rho(:,:);
    rho(ii,1) = rho(ii,2);                                       % Using neighbor point's density
    u(ii,1,:) = 0;
    rho(ii,Ny) = rho(ii,Ny-1);                                       % Using neighbor point's density
    u(ii,Ny,:) = 0;
    rho(1,jj) = rho(2,jj);                                     
    u(1,jj,:) = 0;
     rho(Nx,jj) = rho(Nx-1,jj);                                     
    u(Nx,jj,:) = 0;
    feq_w=zeros(1,Ny-1);feq_f=zeros(1,Ny-1);
 %}
 %上下
 %
flag=1;
    edotu=ev(1,flag)*u(:,Ny-1,1)+ev(2,flag)*u(:,Ny-1,2);usq=u(:,Ny-1,1).^2+u(:,Ny-1,2).^2;
    feq_f=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,Ny,flag)=feq_w+(f(:,Ny-1,flag)-feq_f);
    edotu=ev(1,flag)*u(:,2,1)+ev(2,flag)*u(:,2,2);usq=u(:,2,1).^2+u(:,2,2).^2;
    feq_f=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,1,flag)=feq_w+(f(:,2,flag)-feq_f);
flag=2;
    edotu=ev(1,flag)*u(:,Ny-1,1)+ev(2,flag)*u(:,Ny-1,2);usq=u(:,Ny-1,1).^2+u(:,Ny-1,2).^2;
    feq_f=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,Ny,flag)=feq_w+(f(:,Ny-1,flag)-feq_f);
    edotu=ev(1,flag)*u(:,2,1)+ev(2,flag)*u(:,2,2);usq=u(:,2,1).^2+u(:,2,2).^2;
    feq_f=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,1,flag)=feq_w+(f(:,2,flag)-feq_f);
flag=3;
    edotu=ev(1,flag)*u(:,Ny-1,1)+ev(2,flag)*u(:,Ny-1,2);usq=u(:,Ny-1,1).^2+u(:,Ny-1,2).^2;
    feq_f=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,Ny,flag)=feq_w+(f(:,Ny-1,flag)-feq_f);
    edotu=ev(1,flag)*u(:,2,1)+ev(2,flag)*u(:,2,2);usq=u(:,2,1).^2+u(:,2,2).^2;
    feq_f=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,1,flag)=feq_w+(f(:,2,flag)-feq_f);
flag=4;
    edotu=ev(1,flag)*u(:,Ny-1,1)+ev(2,flag)*u(:,Ny-1,2);usq=u(:,Ny-1,1).^2+u(:,Ny-1,2).^2;
    feq_f=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,Ny,flag)=feq_w+(f(:,Ny-1,flag)-feq_f);
    edotu=ev(1,flag)*u(:,2,1)+ev(2,flag)*u(:,2,2);usq=u(:,2,1).^2+u(:,2,2).^2;
    feq_f=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,1,flag)=feq_w+(f(:,2,flag)-feq_f);
flag=5;
    edotu=ev(1,flag)*u(:,Ny-1,1)+ev(2,flag)*u(:,Ny-1,2);usq=u(:,Ny-1,1).^2+u(:,Ny-1,2).^2;
    feq_f=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,Ny,flag)=feq_w+(f(:,Ny-1,flag)-feq_f);
    edotu=ev(1,flag)*u(:,2,1)+ev(2,flag)*u(:,2,2);usq=u(:,2,1).^2+u(:,2,2).^2;
    feq_f=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,1,flag)=feq_w+(f(:,2,flag)-feq_f);
flag=6;
    edotu=ev(1,flag)*u(:,Ny-1,1)+ev(2,flag)*u(:,Ny-1,2);usq=u(:,Ny-1,1).^2+u(:,Ny-1,2).^2;
    feq_f=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,Ny,flag)=feq_w+(f(:,Ny-1,flag)-feq_f);
    edotu=ev(1,flag)*u(:,2,1)+ev(2,flag)*u(:,2,2);usq=u(:,2,1).^2+u(:,2,2).^2;
    feq_f=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,1,flag)=feq_w+(f(:,2,flag)-feq_f);
flag=7;
    edotu=ev(1,flag)*u(:,Ny-1,1)+ev(2,flag)*u(:,Ny-1,2);usq=u(:,Ny-1,1).^2+u(:,Ny-1,2).^2;
    feq_f=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,Ny,flag)=feq_w+(f(:,Ny-1,flag)-feq_f);
    edotu=ev(1,flag)*u(:,2,1)+ev(2,flag)*u(:,2,2);usq=u(:,2,1).^2+u(:,2,2).^2;
    feq_f=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,1,flag)=feq_w+(f(:,2,flag)-feq_f);
flag=8;
    edotu=ev(1,flag)*u(:,Ny-1,1)+ev(2,flag)*u(:,Ny-1,2);usq=u(:,Ny-1,1).^2+u(:,Ny-1,2).^2;
    feq_f=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,Ny,flag)=feq_w+(f(:,Ny-1,flag)-feq_f);
    edotu=ev(1,flag)*u(:,2,1)+ev(2,flag)*u(:,2,2);usq=u(:,2,1).^2+u(:,2,2).^2;
    feq_f=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,1,flag)=feq_w+(f(:,2,flag)-feq_f);
flag=9;
    edotu=ev(1,flag)*u(:,Ny-1,1)+ev(2,flag)*u(:,Ny-1,2);usq=u(:,Ny-1,1).^2+u(:,Ny-1,2).^2;
    feq_f=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,Ny-1)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,Ny,flag)=feq_w+(f(:,Ny-1,flag)-feq_f);
    edotu=ev(1,flag)*u(:,2,1)+ev(2,flag)*u(:,2,2);usq=u(:,2,1).^2+u(:,2,2).^2;
    feq_f=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(:,Ny-1,1)));usq=zeros(size(u(:,Ny-1,1))).^2;
    feq_w=rho(:,2)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(:,1,flag)=feq_w+(f(:,2,flag)-feq_f);
%}
%左右
    flag=1;
    edotu=ev(1,flag)*u(Nx-1,:,1)+ev(2,flag)*u(Nx-1,:,2);usq=u(Nx-1,:,1).^2+u(Nx-1,:,2).^2;
    feq_f=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(Nx,:,flag)=feq_w+(f(Nx-1,:,flag)-feq_f);
    edotu=ev(1,flag)*u(2,:,1)+ev(2,flag)*u(2,:,2);usq=u(2,:,1).^2+u(2,:,2).^2;
    feq_f=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(1,:,flag)=feq_w+(f(2,:,flag)-feq_f);
    flag=2;
    edotu=ev(1,flag)*u(Nx-1,:,1)+ev(2,flag)*u(Nx-1,:,2);usq=u(Nx-1,:,1).^2+u(Nx-1,:,2).^2;
    feq_f=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(Nx,:,flag)=feq_w+(f(Nx-1,:,flag)-feq_f);
    edotu=ev(1,flag)*u(2,:,1)+ev(2,flag)*u(2,:,2);usq=u(2,:,1).^2+u(2,:,2).^2;
    feq_f=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(1,:,flag)=feq_w+(f(2,:,flag)-feq_f);
    flag=3;
    edotu=ev(1,flag)*u(Nx-1,:,1)+ev(2,flag)*u(Nx-1,:,2);usq=u(Nx-1,:,1).^2+u(Nx-1,:,2).^2;
    feq_f=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(Nx,:,flag)=feq_w+(f(Nx-1,:,flag)-feq_f);
    edotu=ev(1,flag)*u(2,:,1)+ev(2,flag)*u(2,:,2);usq=u(2,:,1).^2+u(2,:,2).^2;
    feq_f=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(1,:,flag)=feq_w+(f(2,:,flag)-feq_f);
    flag=4;
    edotu=ev(1,flag)*u(Nx-1,:,1)+ev(2,flag)*u(Nx-1,:,2);usq=u(Nx-1,:,1).^2+u(Nx-1,:,2).^2;
    feq_f=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(Nx,:,flag)=feq_w+(f(Nx-1,:,flag)-feq_f);
    edotu=ev(1,flag)*u(2,:,1)+ev(2,flag)*u(2,:,2);usq=u(2,:,1).^2+u(2,:,2).^2;
    feq_f=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(1,:,flag)=feq_w+(f(2,:,flag)-feq_f);
    flag=5;
    edotu=ev(1,flag)*u(Nx-1,:,1)+ev(2,flag)*u(Nx-1,:,2);usq=u(Nx-1,:,1).^2+u(Nx-1,:,2).^2;
    feq_f=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(Nx,:,flag)=feq_w+(f(Nx-1,:,flag)-feq_f);
    edotu=ev(1,flag)*u(2,:,1)+ev(2,flag)*u(2,:,2);usq=u(2,:,1).^2+u(2,:,2).^2;
    feq_f=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(1,:,flag)=feq_w+(f(2,:,flag)-feq_f);
    flag=6;
    edotu=ev(1,flag)*u(Nx-1,:,1)+ev(2,flag)*u(Nx-1,:,2);usq=u(Nx-1,:,1).^2+u(Nx-1,:,2).^2;
    feq_f=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(Nx,:,flag)=feq_w+(f(Nx-1,:,flag)-feq_f);
    edotu=ev(1,flag)*u(2,:,1)+ev(2,flag)*u(2,:,2);usq=u(2,:,1).^2+u(2,:,2).^2;
    feq_f=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(1,:,flag)=feq_w+(f(2,:,flag)-feq_f);
    flag=7;
    edotu=ev(1,flag)*u(Nx-1,:,1)+ev(2,flag)*u(Nx-1,:,2);usq=u(Nx-1,:,1).^2+u(Nx-1,:,2).^2;
    feq_f=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(Nx,:,flag)=feq_w+(f(Nx-1,:,flag)-feq_f);
    edotu=ev(1,flag)*u(2,:,1)+ev(2,flag)*u(2,:,2);usq=u(2,:,1).^2+u(2,:,2).^2;
    feq_f=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(1,:,flag)=feq_w+(f(2,:,flag)-feq_f);
    flag=8;
    edotu=ev(1,flag)*u(Nx-1,:,1)+ev(2,flag)*u(Nx-1,:,2);usq=u(Nx-1,:,1).^2+u(Nx-1,:,2).^2;
    feq_f=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(Nx,:,flag)=feq_w+(f(Nx-1,:,flag)-feq_f);
    edotu=ev(1,flag)*u(2,:,1)+ev(2,flag)*u(2,:,2);usq=u(2,:,1).^2+u(2,:,2).^2;
    feq_f=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(1,:,flag)=feq_w+(f(2,:,flag)-feq_f);
    flag=9;
    edotu=ev(1,flag)*u(Nx-1,:,1)+ev(2,flag)*u(Nx-1,:,2);usq=u(Nx-1,:,1).^2+u(Nx-1,:,2).^2;
    feq_f=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(Nx-1,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(Nx,:,flag)=feq_w+(f(Nx-1,:,flag)-feq_f);
    edotu=ev(1,flag)*u(2,:,1)+ev(2,flag)*u(2,:,2);usq=u(2,:,1).^2+u(2,:,2).^2;
    feq_f=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    edotu=ev(1,flag)*zeros(size(u(Nx-1,:,1)));usq=zeros(size(u(Nx-1,:,1)));
    feq_w=rho(2,:)*omega(flag).*(1+a1*edotu+a2*edotu.^2-a3*usq);
    f(1,:,flag)=feq_w+(f(2,:,flag)-feq_f);

  
    
     %{
     %没有必要？
    %固体内部f1=1，其余为0.
k=1;f_trans=f(:,:,k);f_trans(Flag_cir==0)=omega(k);f(:,:,k)=f_trans;
k=2;f_trans=f(:,:,k);f_trans(Flag_cir==0)=omega(k);f(:,:,k)=f_trans;
k=3;f_trans=f(:,:,k);f_trans(Flag_cir==0)=omega(k);f(:,:,k)=f_trans;
k=4;f_trans=f(:,:,k);f_trans(Flag_cir==0)=omega(k);f(:,:,k)=f_trans;
k=5;f_trans=f(:,:,k);f_trans(Flag_cir==0)=omega(k);f(:,:,k)=f_trans;
k=6;f_trans=f(:,:,k);f_trans(Flag_cir==0)=omega(k);f(:,:,k)=f_trans;
k=7;f_trans=f(:,:,k);f_trans(Flag_cir==0)=omega(k);f(:,:,k)=f_trans;
k=8;f_trans=f(:,:,k);f_trans(Flag_cir==0)=omega(k);f(:,:,k)=f_trans;
k=9;f_trans=f(:,:,k);f_trans(Flag_cir==0)=omega(k);f(:,:,k)=f_trans;
    %}
%5，求宏观u，rho
rho(:,:)=sum(f,3);
u(:,:,1)=f(:,:,1)*ev(1,1)+f(:,:,2)*ev(1,2)+f(:,:,3)*ev(1,3)+f(:,:,4)*ev(1,4)+f(:,:,5)*ev(1,5)+...
        f(:,:,6)*ev(1,6)+f(:,:,7)*ev(1,7)+f(:,:,8)*ev(1,8)+f(:,:,9)*ev(1,9)+F(:,:,1)*delta_t/2;
u(:,:,2)=f(:,:,1)*ev(2,1)+f(:,:,2)*ev(2,2)+f(:,:,3)*ev(2,3)+f(:,:,4)*ev(2,4)+f(:,:,5)*ev(2,5)+...
        f(:,:,6)*ev(2,6)+f(:,:,7)*ev(2,7)+f(:,:,8)*ev(2,8)+f(:,:,9)*ev(2,9)+(F(:,:,2)+G)*delta_t/2;
u(:,:,1)=u(:,:,1)./rho(:,:);
u(:,:,2)=u(:,:,2)./rho(:,:);
%

%传热

uTx=u(:,:,1);uTy=u(:,:,2);
uTx(Flag_squ==1)=0;uTy(Flag_squ==1)=0;
%T=T.*Trc;

f1T=permute(h,[3,1,2]);
mT=MT*f1T(1:5,:);
     mT_eq(1,:,:)=T;
     mT_eq(2,:,:)=T.*uTx;
     mT_eq(3,:,:)=T.*uTy;
     mT_eq(4,:,:)=2*T/3;
     %mT_eq(4,:,:)=1*T/3;
     mT_eq(5,:,:)=0;
mT_temp0=Minv_sT*(mT-mT_eq(1:5,:));
%mS_temp0=MinvS*mS_temp0;
mT_temp=reshape(mT_temp0,5,Nx,Ny);f1T=f1T-mT_temp;    
h=permute(f1T,[2,3,1]);

k=1;h(:,:,k) =circshift(h(:,:,k),[evS(1,k),evS(2,k),0]);
k=2;h(:,:,k) =circshift(h(:,:,k),[evS(1,k),evS(2,k),0]);
k=3;h(:,:,k) =circshift(h(:,:,k),[evS(1,k),evS(2,k),0]);
k=4;h(:,:,k) =circshift(h(:,:,k),[evS(1,k),evS(2,k),0]);
k=5;h(:,:,k) =circshift(h(:,:,k),[evS(1,k),evS(2,k),0]);

%{
%S=sum(s,3);
pcolor(X(:,1:end-1),Y(:,1:end-1),S(1:end-1,:)')
%pcolor(X(:,1:end-1),Y(:,1:end-1),T(1:end-1,:,1)')
    colormap jet;
    colorbar('location','eastoutside')
    shading interp
    pause(2)
%}


%6，s,h边界
h(:,Ny,5)=-1/6-h(:,Ny,3);h(:,1,3)=1/6-h(:,1,5);
h(1,:,2)=h(1,:,4); h(Nx,:,4)=h(Nx,:,2); 
%h(:,Ny,5)=0-h(:,Ny,3);h(:,1,3)=1/3-h(:,1,5);
%h(Nx,:,4)=0-h(Nx,:,2);h(1,:,2)=0-h(1,:,4);
%
T=sum(h,3);
%right,left
    %绝热

%S=sum(s,3)+Sf*delta_t/2;
%S=sum(s,3);
%S(S>1)=1;S(S<0)=0;   %zxf
%S(:,1)=0;S(:,end)=1;
%T(:,1)=1;T(:,end)=0;
T(:,1)=0.5;T(:,end)=-0.5;
%迭代浸没
%{
error=nerror;
while (error>0)
    
    %内边界 %dS/dn=0;
    %修正温差
    dTi=2*((delt2o-delt2i)*S(:))/delta_x^2;
    %内部分散
    linshi=reshape(del2i*dTi*(ds*(R-delta_x)/R),Nx,Ny);
    linshi(Flag_cirin==1)=0;
    Terror(:,:,error)=linshi;
    %3，更新S
    S(:,:)=S(:,:)+lamta*Terror(:,:,error)*delta_t/2;
    error=error-1;
end
Sf=sum(Terror,3)+Sf;
%?
%共轭传热
%修正温差
    dTi=2*((delt2o+kfs*delt2i-(1+kfs)*delt)*T(:))/(1+kfs);
    %dTi=-((delt2o-delt2i)*T(:))/(1+0.5)*mu/Pr;
    
    %内部分散
    linshi=reshape(del*dTi*ds,Nx,Ny);
    %Terror(:,:,error)=Terror(:,:,error)+(error==5)*linshi;
     %
    %3，更新T
    %T(:,:)=T(:,:)+lamta*Terror(:,:,error)*delta_t/2;
     T(:,:)=T(:,:)+2*linshi*delta_t/2;
%}
%S(S>0.5)=0.5;S(S<-0.5)=-0.5;
T(T>0.5)=0.5;T(T<-0.5)=-0.5;
%T(T>1)=1;T(T<0)=-0;
%}
%7.5,示踪粒子

%{
ux=u(:,:,1);ux=ux(:);uy=u(:,:,2);uy=uy(:);
szx=round(shizongx(:,t)/0.005)+1;szy=round(shizongy(:,t)/0.005)+1;
shizongu=ux((szy-1)*(Nx)+szx);shizongv=uy((szy-1)*(Nx)+szx);
shizongx(:,t+1)=shizongu*delta_t+shizongx(:,t);
shizongy(:,t+1)=shizongv*delta_t+shizongy(:,t);
szf=shizongx(:,t+1);szf(szf>2)=szf(szf>2)-2;szf(szf<0)=szf(szf<0)+2;shizongx(:,t+1)=szf;
szf=shizongy(:,t+1);szf(szf>2)=2;szf(szf<0)=0;shizongy(:,t+1)=szf;
%}

%8，作图，保存数据
%
if( mod(t,PlotInterval) ==0 )
    clf
    %pS=S;pS(S>1)=1+(S(S>1)-1)*100;   %zxf
    pcolor(X(:,1:end-1),Y(:,1:end-1),T(1:end-1,:)')
 
    colormap jet;
    %
    %{
    pcolor(X,Y,(u(1:Nx,:,1).^2+u(1:Nx,:,2).^2)')
    %colormap(mycolor);
     streamslice(X,Y,u(:,:,1)',u(:,:,2)')
    %}
    colorbar('location','eastoutside')
    shading interp
    axis equal
    axis([0,Lx,0,Ly])
    hold on
   
%    fill(Xlag2,Ylag2,'k');
    box on
    %{
    %plot(xPos,yPos,'o','MarkerFaceColor','k','MarkerSize',400*0.046)
    %cir_x=repmat(xPos,629,1)+repmat(0.046*cos((0:0.01:2*pi)'),1,45);
    %cir_y=repmat(yPos,629,1)+repmat(0.046*sin((0:0.01:2*pi)'),1,45);
    %}
    %plot([sxb(:),sxe(:),sxe(:),sxb(:),sxb(:)]',[syb(:),syb(:),sye(:),sye(:),syb(:)]','k')
    %fill(cir_x,cir_y,'w');
    fill([sxb(:),sxe(:),sxe(:),sxb(:),sxb(:)]',[syb(:),syb(:),sye(:),sye(:),syb(:)]','k')
    hold off
    %title(num2str(t))
    pause(0.01)
    %frame = getframe;%2
    %writeVideo(writerObj,frame);%2
end

if (mod(t,10) ==0)
Re(t/10,1)= (sum(sum((uTx.^2+uTy.^2)))/Area).^0.5/mu;
Nu(t/10,1)=abs(-sum(-.5-T(:,end-1))/Lx)*(Nx-1)/Nx;
end


if( mod(t,OutputsInterval) ==0)
       %file_name = FieldOutputs(rho, u, t, files_name);
       %fprintf('Time Step = %d; data saved as: %s \n', t,file_name);
       save([num2str(t),'RB.mat'],'t','u','T','rho');   
       if(mod(t,50000) ==0)
           save([num2str(t),'main.mat']); 
       end  
end

%


end

%}

%close(writerObj);%3
fprintf('END\n');
%save t=20000.mat
%save Ra=5700.mat