%%
% Nume si prenume: Cojocaru Darius-Andrei
%

clearvars
clc

%% Magic numbers (replace with received numbers)
m = 6; 
n = 8; 

%% Process data (fixed, do not modify)
c1 = (1000+n*300)/10000;
c2 = (1.15+2*(m+n/10)/20);
a1 = 2*c2*c1;
a2 = c1;
b0 = (1.2+m+n)/5.5;

rng(m+10*n)
x0_slx = [2*(m/2+rand(1)*m/5); m*(n/20+rand(1)*n/100)];

%% Experiment setup (fixed, do not modify)
Ts = 10*c2/c1/1e4*1.5; % fundamental step size
Tfin = 30*c2/c1*10; % simulation duration

gain = 10;
umin = 0; umax = gain; % input saturation
ymin = 0; ymax = b0*gain/1.5; % output saturation

whtn_pow_in = 1e-6*5*(((m-1)*8+n/2)/5)/2*6/8; % input white noise power and sampling time
whtn_Ts_in = Ts*3;
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(10); % input quantizer (DAC)

whtn_pow_out = 1e-5*5*(((m-1)*25+n/2)/5)*6/80*(0.5+0.3*(m-2)); % output white noise power and sampling time
whtn_Ts_out = Ts*5;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(9); % output quantizer (ADC)

u_op_region = (m/2+n/5)/2; % operating point

%% Input setup (can be changed/replaced/deleted)
wf=1/10.03 %%1/T1, T1 dominant, usor de aproximat din raspunsul la treapta
fmin=wf/2/pi/10
fmax=wf/2/pi*5 %scade-l
Ain=1.48 %arbitrar, suficient de mare incat sa depaseasca zgomotul

%% Data acquisition (use t, u, y to perform system identification)
out = sim("circuit_hidraulic_R2022b.slx");

t = out.tout;
u = out.u;
y = out.y;

u=medfilt1(u,30);
y=medfilt1(y,30);

plot(t,u,t,y)
shg

%% System identification
y_st=(6.6572+6.00947)/2
u_st=(3.75+.874023)/2
K=y_st/u_st
w1=pi/(319.066-305.018)
DeltaT1=312.284-305.018
phi1=rad2deg(-w1*DeltaT1)%% sa iasa 90 de grade
Ay=(7.91667-4.76799)/2
Au=1.48
M=Ay/Au
Im=-M 
wn=w1
zeta=-K/2/Im 
H=tf(K*wn^2,[1,2*zeta*wn,wn^2])
zpk(H)
T1=1/.4693
T2=1/.1066

A=[0,1;-1/T1/T2,-(1/T1+1/T2)];
B=[0;K/T1/T2];
C=[1,0];
D=0;
sys=ss(A,B,C,D)
ysim2=lsim(sys,u,t,[11,.2]);

figure
plot(t,u,t,y,t,ysim2)

J=1/sqrt(length(t))*norm(y-ysim2)
eMPN=norm(y-ysim2)/norm(y-mean(y))*100