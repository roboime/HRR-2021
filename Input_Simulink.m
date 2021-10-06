  %Humanide Input
%% Contact/friction parameters
contact_stiffness = 25000;
contact_damping = 100;
mu_k = 0.6;
mu_s = 0.8;
mu_vth = 0.1;
height_plane = 0.025;
dz=-0.333;
plane_x = 5;
plane_y = 1;
contact_point_radius = 1e-4;
%% Outros
tq_max_mt=40;
Tmax=40;
g=9.80665;
zzmp=0.012446-dz;
chao=[plane_x plane_y height_plane];
damp_f=5e-4;
filter=0.001;
%t_e=13*t0;
%ta=7*t0-0.2;
%massabat=0.045;
%Ibat=[3.55875  10.59 11.71875];
%Pbat=[0 0 0];
%% PID
P=300;
I=5;
D=20;
RT=5e4;
%P=39.42;
%I=3.077;
%D=0.77;
%% Geometria
l1=67;
l2=57;
l3=53.3038;