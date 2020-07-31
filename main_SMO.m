%-----------------------------------------------------------------------
%   This is the code for the following paper:
%
%   Q. He, W. Zhang, P. Lu and J. Liu
%   Performance comparison of representative model-based fault 
%   reconstruction algorithms for aircraft sensor fault detection 
%   and diagnosis.
%
%   if you find this code helpful, please cite the above paper.
%   
%   If you have more questions, please contact the correponding author:
%   Peng Lu   
%   Email:  peng.lu@polyu.edu.hk
%   Released in 2019
%----------------------------------------------------------------------

%% START

clc
clear
close all

tic;

gk = 9.80665;
r2d = 180/pi;

GPS_rate = [1; 10; 100];
num = 3;
S_hz = GPS_rate(num);
load new_sim_Ma_06_H_1000_notur.mat;

%% organize data

t_s = 1;
t_e = size(Out_nonoise,1);

% true states
Vtas   =  Out_nonoise(t_s:t_e,1)';
alpha  =  Out_nonoise(t_s:t_e,2)';
beta   =  Out_nonoise(t_s:t_e,3)';
pb     =  Out_nonoise(t_s:t_e,4)';
qb     =  Out_nonoise(t_s:t_e,5)';
rb     =  Out_nonoise(t_s:t_e,6)';
phi    =  Out_nonoise(t_s:t_e,9)';
theta  =  Out_nonoise(t_s:t_e,8)';
psi    =  Out_nonoise(t_s:t_e,7)';
xe     =  Out_nonoise(t_s:t_e,10)';
ye     =  Out_nonoise(t_s:t_e,11)';
ze     =  Out_nonoise(t_s:t_e,12)'; % opposite sign of he
u_b    =  Out_nonoise(t_s:t_e,13)';
v_b    =  Out_nonoise(t_s:t_e,14)';
w_b    =  Out_nonoise(t_s:t_e,15)';
u_n    =  Out_nonoise(t_s:t_e,16)';
v_n    =  Out_nonoise(t_s:t_e,17)';
w_n    =  Out_nonoise(t_s:t_e,18)'; % the opposite sign of hdot

% A_x    =  Out_nonoise(t_s:t_e,32)';
% A_y    =  Out_nonoise(t_s:t_e,33)';
% A_z    =  Out_nonoise(t_s:t_e,34)';

% % without noise
% A_x    =  Ax';
% A_y    =  Ay';
% A_z    =  Az';

% with noise
A_x    =  squeeze(Ax_m)';
A_y    =  squeeze(Ay_m)';
A_z    =  squeeze(Az_m)';

% turbulence
ug_b     = Tur_b(t_s:t_e,1)';
vg_b     = Tur_b(t_s:t_e,2)';
wg_b     = Tur_b(t_s:t_e,3)';
% p_tur    = Tur_b(t_s:t_e,4)';
 
% derivative of the turbulence
ug_b_dot     = Tur_b_dot(t_s:t_e,1)';
vg_b_dot     = Tur_b_dot(t_s:t_e,2)';
wg_b_dot     = Tur_b_dot(t_s:t_e,3)';

% NOTE the derivative comes from last step!!! 
% otherwise will cause a delay
ug_b_dot      =  [0 ug_b_dot(1:length(ug_b_dot)-1)];
vg_b_dot      =  [0 vg_b_dot(1:length(ug_b_dot)-1)];
wg_b_dot      =  [0 wg_b_dot(1:length(ug_b_dot)-1)];

% grond velocity in the body axis 
u_n_b = u_n.*cos(theta).*cos(psi)+v_n.*cos(theta).*sin(psi)-w_n.*sin(theta);
v_n_b = u_n.*(sin(theta).*cos(psi).*sin(phi)-sin(psi).*cos(phi))...
    +v_n.*(sin(theta).*sin(psi).*sin(phi)+cos(psi).*cos(phi))...
    +w_n.*cos(theta).*sin(phi);
w_n_b = u_n.*(sin(theta).*cos(psi).*cos(phi)+sin(psi).*sin(phi))...
    +v_n.*(sin(theta).*sin(psi).*cos(phi)-cos(psi).*sin(phi))...
    +w_n.*cos(theta).*cos(phi);

% lower the dimension
Out_noise=squeeze(Out_noise)';

% measurenments
p_m         =  Out_noise(t_s:t_e,4)';
q_m         =  Out_noise(t_s:t_e,5)';
r_m         =  Out_noise(t_s:t_e,6)';
Ax_m        =  squeeze(Ax_m)'; % due to the simulation
Ay_m        =  squeeze(Ay_m)';
Az_m        =  squeeze(Az_m)';
Ax_m        =  Ax_m(t_s:t_e); % due to the simulation
Ay_m        =  Ay_m(t_s:t_e);
Az_m        =  Az_m(t_s:t_e);
Vt_m        =  Out_noise(t_s:t_e,1)';
alpha_m     =  Out_noise(t_s:t_e,2)';
beta_m      =  Out_noise(t_s:t_e,3)';
phi_m       =  Out_noise(t_s:t_e,9)';
theta_m     =  Out_noise(t_s:t_e,8)';
psi_m       =  Out_noise(t_s:t_e,7)';
ze_m        =  Out_noise(t_s:t_e,12)'; % opposite sign of he_m
% ze          =  Out_noise(t_s:t_e,17)';
xe_m        =  Out_noise(t_s:t_e,10)';
ye_m        =  Out_noise(t_s:t_e,11)';
% u_b_m       =  Out_noise(t_s:t_e,13)'; % not realistic
% v_b_m       =  Out_noise(t_s:t_e,14)'; % not realistic
% w_b_m       =  Out_noise(t_s:t_e,15)'; % not realistic
u_n_m       =  Out_noise(t_s:t_e,16)';
v_n_m       =  Out_noise(t_s:t_e,17)';
w_n_m       =  Out_noise(t_s:t_e,18)'; % Vz= -hedot, Vz and wn are in the same direction

%-------------------------- different sampling rate -------------------
u_n_m0 = u_n_m;
v_n_m0 = v_n_m;
w_n_m0 = w_n_m;
dtu_g = 1/S_hz;
dtv_g = 1/S_hz;
dtw_g = 1/S_hz;
u_n_m = sampling_g(u_n_m0,d_t,dtu_g);
v_n_m = sampling_g(v_n_m0,d_t,dtv_g);
w_n_m = sampling_g(w_n_m0,d_t,dtw_g);

%-------------------------- different sampling rate -------------------
UN_m = u_n_m;
UE_m = v_n_m;
UD_m = w_n_m;

%% add noise to the true state

while(true)
    
    noise_flag =input('Please select the noise\n1.standard noise\n2.amplified noise\n');
    
    if noise_flag==1
        a=1e-4;b=3e-8;
        break;
    end
    
    if noise_flag==2
        a=1e-3;b=3e-7;
        break;
    end
end

Ax_m = Ax'+sqrt(a)*randn(1,10001);
Ay_m = Ay'+sqrt(a)*randn(1,10001);
Az_m = Az'+sqrt(a)*randn(1,10001);
p_m = pb+sqrt(b)*randn(1,10001);
q_m = qb+sqrt(b)*randn(1,10001);
r_m = rb+sqrt(b)*randn(1,10001);
UN_m = u_n +sqrt(a)*randn(1,10001);
UE_m = v_n +sqrt(a)*randn(1,10001);
UD_m = w_n +sqrt(a)*randn(1,10001);
phi_m = phi +sqrt(b)*randn(1,10001);
theta_m = theta +sqrt(b)*randn(1,10001);
psi_m = psi +sqrt(b)*randn(1,10001);

%% add IMU faults  

u  =  [Ax_m; Ay_m; Az_m; p_m; q_m; r_m];
[dim_u,gen] = size(u);
lameda = [ 0; 0; 0; 0; 0; 0 ]*ones(1,gen);
x_real = [ u_n; v_n; w_n; phi; theta; psi;];

u_origin = u;
faults_in = zeros(6,gen);
fy = 0.001;

faults_in(1,1001:3000) = faults_in(1,1001:3000) + 1*2;
faults_in(2,1001:3000) = faults_in(2,1001:3000) + fy*(1:2000);
faults_in(3,1001:3000) = faults_in(3,1001:3000) + sin(0.005*pi*(1001:3000))*2;
faults_in(4,1001:3000) = faults_in(4,1001:3000) + 1/57.3*2;
faults_in(5,1001:3000) = faults_in(5,1001:3000) + 1/57.3*fy*(1:2000);
faults_in(6,1001:3000) = faults_in(6,1001:3000) + 1/57.3*sin(0.005*pi*(1:2000))*2;

u(1,:) = u(1,:) + faults_in(1,:);
u(2,:) = u(2,:) + faults_in(2,:);
u(3,:) = u(3,:) + faults_in(3,:);
u(4,:) = u(4,:) + faults_in(4,:);
u(5,:) = u(5,:) + faults_in(5,:);
u(6,:) = u(6,:) + faults_in(6,:);

%% Sliding mode observer

z_real  =  [ UN_m; UE_m; UD_m; phi_m; theta_m; psi_m;];
x_ob_0 = [ UN_m(2); UE_m(2); UD_m(2); phi_m(1); theta_m(1); psi_m(1); ];
[dim_sys,~] = size(x_ob_0);
[dim_out,~] = size(z_real);
x_ob = zeros(dim_sys,gen);
z_ob = zeros(dim_out,gen);
x_ob_dot = zeros(dim_sys,gen); 
v_eq = zeros(dim_sys,gen);
ey = zeros(dim_out,gen);
f_ob = zeros(dim_u,gen);
x_ob_filter = zeros(dim_sys,gen);
z_ob_filter = zeros(dim_out,gen);
residual = zeros(dim_out,gen);
norm_residual = zeros(1,gen);
err_x = zeros(dim_sys,gen);
norm_error_x = zeros(1,gen);
R = zeros(dim_sys,gen);
Is = zeros(dim_sys,gen);
f_ob_filter_dot = zeros(dim_u,gen);
f_ob_filter = zeros(dim_u,gen);
dim_fi = dim_u;
n = dim_sys;
p = dim_out;
m = dim_u;
q = dim_u;

A = F_lin(x_real(:,1),u(:,1),9.80665);
B = - G_lin(x_real(:,1));
E = - E_lin(x_real(:,1));
C = eye(dim_out);
G = B; 

eig(A);
T_O = ctrb(A,C);
rank(T_O);

mu_i = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1];
mu_i_one = 0.1;
K = 10; % should be scaled together with mu_i
% design L
Pole = [-3.01, -3.02, -3.03, -3.04, -3.05, -3.06];
% Pole = [-4, -5, -6, -7, -8, -9];
% Pole = [-11, -12, -13, -14, -15, -16];
[L, pre, message] = place(A, C, Pole);
A_new = A - L*C;
tau = 0.15; % low pass filter
f_ob_filter_0 = 0;

for k = 1:5000

%------------- use time-varying matrices
%     A = F_lin(x_ob_0,u(:,k),9.80665);
%     B = - G_lin(x_ob_0);
%     E = - E_lin(x_ob_0);
%     C = eye(dim_out);

    z_ob(:,k) = C * x_ob_0;
    inv_E = eye(dim_fi)/E;
    err_x(:,k) = ( eye(dim_out)/C ) * z_real(:,k) - x_ob_0;
    R(:,k) = inv_E * err_x(:,k);
    Is(:,k) = tanh( (R(:,k)+ 0.0) ./ mu_i);
    % observer
    F=fx(x_ob_0,u(:,k),gk);    
    x_ob_dot(:,k) = (A - L*C) * x_ob_0 + L* z_real(:,k) + B * u(:,k) + K * Is(:,k);
    x_ob_dot(3,k) = x_ob_dot(3,k) + gk; % include the g
    x_ob(:,k) = x_ob_0 + x_ob_dot(:,k)* d_t;
    x_ob_0 = x_ob(:,k);    
    % fault estimation
    f_ob(:,k) = inv_E * (K * Is(:,k));
    f_ob(:,k) = - f_ob(:,k);
    % low pass filter
    f_ob_filter_dot(:,k) = 1/tau* ( f_ob(:,k) - f_ob_filter_0 );
    f_ob_filter(:,k) = f_ob_filter_0 + f_ob_filter_dot(:,k) * d_t;
    f_ob_filter_0 = f_ob_filter(:,k);
          
end
toc

%% plot
Time=d_t*(1:5000);

figure;
subplot(311);hold on;plot(Time,f_ob_filter(1,1:k),'b','linewidth',1.5);plot(Time,faults_in(1,1:k),'r--','linewidth',1.5);
ylabel('f_{Ax} (m/s^2)','fontsize',13);grid;
set(gca,'xlim',[0 d_t*k],'fontsize',13); set(gca,'ylim',[-0.8 2.5]);
subplot(312);hold on;plot(Time,f_ob_filter(2,1:k),'b','linewidth',1.5);plot(Time,faults_in(2,1:k),'r--','linewidth',1.5); 
ylabel('f_{Ay} (m/s^2)','fontsize',13);grid;
set(gca,'xlim',[0 d_t*k],'fontsize',13); set(gca,'ylim',[-0.8 2.2]);
subplot(313);hold on;plot(Time,f_ob_filter(3,1:k),'b','linewidth',1.5);plot(Time,faults_in(3,1:k),'r--','linewidth',1.5);
ylabel('f_{Az} (m/s^2)','fontsize',13);grid;
set(gca,'xlim',[0 d_t*k],'fontsize',13); set(gca,'ylim',[-2.2 2.2]);
legend('estimation','true');
title('time (s)','fontsize',13)   
    
figure;
subplot(311);hold on;plot(Time,f_ob_filter(4,1:k),'b','linewidth',1.5);plot(Time,faults_in(4,1:k),'r--','linewidth',1.5);
ylabel('f_{p} (rad/s)','fontsize',13);grid;
set(gca,'xlim',[0 d_t*k],'fontsize',13); set(gca,'ylim',[-0.02 0.05]);
subplot(312);hold on;plot(Time,f_ob_filter(5,1:k),'b','linewidth',1.5);plot(Time,faults_in(5,1:k),'r--','linewidth',1.5);
ylabel('f_{q} (rad/s)','fontsize',13);grid;
set(gca,'xlim',[0 d_t*k],'fontsize',13); set(gca,'ylim',[-0.015 0.04]);
subplot(313);hold on;plot(Time,f_ob_filter(6,1:k),'b','linewidth',1.5);plot(Time,faults_in(6,1:k),'r--','linewidth',1.5);
ylabel('f_{r} (rad/s)','fontsize',13);grid;
set(gca,'xlim',[0 d_t*k],'fontsize',13);  set(gca,'ylim',[-0.04 0.04]);
legend('estimation','true');
title('time (s)','fontsize',13)
%% END