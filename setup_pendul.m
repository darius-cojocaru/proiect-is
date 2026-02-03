%%
% Nume si prenume: Cojocaru Darius-Andrei
%

clearvars
clc

%% Magic numbers (replace with received numbers)
m = 6;
n = 8;

%% Process data and experiment setup (fixed, do not modify)
Ts = 500e-6; % fundamental step size

u_star = 1.2+n*0.075;
delta = 0.125;
delta_spab = 0.075;

umin = -5; umax = 5; % input saturation
ymin = -100; ymax = 100; % output saturation

g = 9.81;
% pendulum parameters
M = 0.8-n/48;
l = 1.2-m/24;
b = 0.3+m/24;
% measurement
c1 = 180/pi;
c2 = 4+n/2;

% (theta0,omega0)
rng(m+10*n)
x0_slx = [(n+3)/50,(-1)^(n+1)*m/20];

% input white noise power and sampling time
whtn_pow_in = 1e-10*(Ts*1e4)/2; 
whtn_Ts_in = Ts*2;
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(13); % input quantizer (DAC)

% output white noise power and sampling time
whtn_pow_out = 1e-3*Ts; 
whtn_Ts_out = Ts*2;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(13); % output quantizer (ADC)

meas_rep = round(7+n/2); % data acquisition hardware sampling limitation

%% Input setup (can be changed/replaced/deleted)
t1=12;%timpul aplanare ci
tr=0.6*2;%timpul primului maxim
N=5;%nr biti registru spab
p=round(tr/N/Ts)%divizor de frecventa spab
DeltaT=p*(2^N-1)*Ts*3;%durata semnalululi
[input_LUT_dSpace,Tfin]=generate_input_signal(Ts,t1,DeltaT,N,p,u_star,delta,delta_spab)
%% Data acquisition (use t, u, y to perform system identification)
out = sim("pendul_R2022b.slx");

t = out.tout;
u = out.u;
y = out.y;

u=medfilt1(u,3);
y=medfilt1(y,3);

subplot(211)
plot(t,u)
subplot(212)
plot(t,y)
shg

%% System identification
i1=45596;%%2 de sus 2 de jos de la trapez
i2=100277;
i3=121271;
i4=176705;

Nr=15;

t_id=t(i1:Nr:i2);
u_id=u(i1:Nr:i2);
u_id=u_id-mean(u_id);
y_id=y(i1:Nr:i2);
y_id=y_id-mean(y_id);

t_vd=t(i3:Nr:i4);
u_vd=u(i3:Nr:i4);
u_vd=u_vd-mean(u_vd);
y_vd=y(i3:Nr:i4);
y_vd=y_vd-mean(y_vd);
figure
subplot(221)
plot(t_id,u_id)
subplot(223)
plot(t_id,y_id)
subplot(222)
plot(t_vd,u_vd)
subplot(224)
plot(t_vd,y_vd)

dat_id=iddata(y_id,u_id,t_id(2)-t_id(1))
dat_vd=iddata(y_vd,u_vd,t_vd(2)-t_vd(1))
%%
% gresit
model_arx=arx(dat_id,[2 2 1])
figure,resid(model_arx,dat_vd)
figure,compare(model_arx,dat_vd)
%%
%e bun fit de 95 trece si testul de autocorelatie 
model   _armax=armax(dat_id,[2 2 5 1]) %na poli, nb zerouri+1, nk tacti intarziere, nc polinomul care mod zgomotul
figure,resid(model_armax,dat_vd)
figure,compare(model_armax,dat_vd)
%%
%e bun fit de 95 trece testul de intercorelatie
model_oe=oe(dat_vd,[2 2 1])%nb nf nk nf=na
figure,resid(model_oe,dat_id)
figure,compare(model_oe,dat_id)
%%
%fit bun 95 trece macar un test
model_ssest=ssest(dat_id,2)%ord sist
figure,resid(model_ssest,dat_vd)
figure,compare(model_ssest,dat_vd)
sys_tf=tf(model_ssest)
