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

d_t = 0.01;
gk = 9.80665;
gen=5000;
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
%--------------  NOTE the derivative comes from last step!!! 
%--------------  otherwise will cause a delay
ug_b_dot      =  [0 ug_b_dot(1:length(ug_b_dot)-1)];
vg_b_dot      =  [0 vg_b_dot(1:length(ug_b_dot)-1)];
wg_b_dot      =  [0 wg_b_dot(1:length(ug_b_dot)-1)];
%-----------------------------------------------------------

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

IMU  =  [Ax_m(1:gen); Ay_m(1:gen); Az_m(1:gen); p_m(1:gen); q_m(1:gen); r_m(1:gen)];
MEA  =  [ UN_m; UE_m; UD_m; phi_m; theta_m; psi_m;];

IMU_origin = IMU;
faults_in = zeros(6,gen);
fy = 0.001;

faults_in(1,1001:3000) = faults_in(1,1001:3000) + 1 * 2;
faults_in(2,1001:3000) = faults_in(2,1001:3000) + fy*(1:2000);
faults_in(3,1001:3000) = faults_in(3,1001:3000) + sin(0.005*pi*(1001:3000))*2;
faults_in(4,1001:3000) = faults_in(4,1001:3000) + 1/57.3*2;
faults_in(5,1001:3000) = faults_in(5,1001:3000) + 1/57.3*fy*(1:2000);
faults_in(6,1001:3000) = faults_in(6,1001:3000) + 1/57.3*sin(0.005*pi*(1:2000))*2;

IMU(1,:) = IMU(1,:) + faults_in(1,:);
IMU(2,:) = IMU(2,:) + faults_in(2,:);
IMU(3,:) = IMU(3,:) + faults_in(3,:);
IMU(4,:) = IMU(4,:) + faults_in(4,:);
IMU(5,:) = IMU(5,:) + faults_in(5,:);
IMU(6,:) = IMU(6,:) + faults_in(6,:);

%% Nonlinear disturbance observer

p = zeros(6,1);
z = zeros(6,gen);
d_est = zeros(6,gen);

tau = 0.15;
f_ob_filter_0 = 0;
f_ob_filter_dot = zeros(6,5000);
f_ob_filter = zeros(6,5000);


 for k=1:gen

     x = MEA(:,k);
     
     [G1, G2] = Matrix_get(x);
     
     A_ob = zeros(6,6);   
     %Pole = [-2, -2, -2, -2, -2, -2];
     %Pole = [-20, -20, -20, -20, -20, -20];
     %Pole = [-10, -10, -10, -10, -10, -10];
     %Pole = [-100, -100, -100, -100, -100, -100];
     Pole = [-4.1,-4.2,-4.3,-4.4,-4.5,-4.6];
     %Pole = [-15,-15,-15,-15,-15,-15];
     L = place(A_ob,G2,Pole);
          
     p(1,:)=L(1,1)*x(1,1)+L(1,2)*x(2,1)+L(1,3)*x(3,1)+L(1,4)*x(4,1)+L(1,5)*x(5,1)+L(1,6)*x(6,1);
     p(2,:)=L(2,1)*x(1,1)+L(2,2)*x(2,1)+L(2,3)*x(3,1)+L(2,4)*x(4,1)+L(2,5)*x(5,1)+L(2,6)*x(6,1);
     p(3,:)=L(3,1)*x(1,1)+L(3,2)*x(2,1)+L(3,3)*x(3,1)+L(3,4)*x(4,1)+L(3,5)*x(5,1)+L(3,6)*x(6,1);
     p(4,:)=L(4,1)*x(1,1)+L(4,2)*x(2,1)+L(4,3)*x(3,1)+L(4,4)*x(4,1)+L(4,5)*x(5,1)+L(4,6)*x(6,1);
     p(5,:)=L(5,1)*x(1,1)+L(5,2)*x(2,1)+L(5,3)*x(3,1)+L(5,4)*x(4,1)+L(5,5)*x(5,1)+L(5,6)*x(6,1);
     p(6,:)=L(6,1)*x(1,1)+L(6,2)*x(2,1)+L(6,3)*x(3,1)+L(6,4)*x(4,1)+L(6,5)*x(5,1)+L(6,6)*x(6,1);
     
     Fx = [0;0;gk;0;0;0];
     
     dot_z = -L*G2*z(:,k)-L*(G2*p+Fx+G1*IMU(:,k));
     
     z(:,k+1) = z(:,k) + dot_z*d_t;
     d_est(:,k) = z(:,k+1) + p;

%      % low pass filter
%      f_ob_filter_dot(:,k) = 1/tau* ( d_est(:,k) - f_ob_filter_0 );
%      f_ob_filter(:,k) = f_ob_filter_0 + f_ob_filter_dot(:,k) * d_t;
%      f_ob_filter_0 = f_ob_filter(:,k);
       
 end
toc;

%% plot

Time=d_t*(1:gen);

figure;
subplot(311);hold on;plot(Time,d_est(1,1:k),'b','linewidth',1.5);plot(Time,faults_in(1,1:k),'r--','linewidth',1.5); 
ylabel('f_{Ax} (m/s^2)','fontsize',13);grid;
set(gca,'xlim',[0 d_t*k],'fontsize',13); set(gca,'ylim',[-1.5 4]);
subplot(312);hold on;plot(Time,d_est(2,1:k),'b','linewidth',1.5); plot(Time,faults_in(2,1:k),'r--','linewidth',1.5); 
ylabel('f_{Ay} (m/s^2)','fontsize',13);grid;
set(gca,'xlim',[0 d_t*k],'fontsize',13); set(gca,'ylim',[-2 3]);%set(gca,'ylim',[-1 2.4]);
subplot(313);hold on;plot(Time,d_est(3,1:k),'b','linewidth',1.5); plot(Time,faults_in(3,1:k),'r--','linewidth',1.5); 
ylabel('f_{Az} (m/s^2)','fontsize',13);grid;
set(gca,'xlim',[0 d_t*k],'fontsize',13); set(gca,'ylim',[-3 3]);
legend('estimation','true');
title('time (s)','fontsize',13)   
    
figure;
subplot(311);hold on;plot(Time,d_est(4,1:k),'b','linewidth',1.5);plot(Time,faults_in(4,1:k),'r--','linewidth',1.5);
ylabel('f_{p} (rad/s)','fontsize',13);grid;
set(gca,'xlim',[0 d_t*k],'fontsize',13);  set(gca,'ylim',[-0.02 0.05]);
subplot(312);hold on;plot(Time,d_est(5,1:k),'b','linewidth',1.5);plot(Time,faults_in(5,1:k),'r--','linewidth',1.5);
ylabel('f_{q} (rad/s)','fontsize',13);grid;
set(gca,'xlim',[0 d_t*k],'fontsize',13);  set(gca,'ylim',[-0.02 0.04]);
subplot(313);hold on;plot(Time,d_est(6,1:k),'b','linewidth',1.5);plot(Time,faults_in(6,1:k),'r--','linewidth',1.5);
ylabel('f_{r} (rad/s)','fontsize',13);grid;
set(gca,'xlim',[0 d_t*k],'fontsize',13);  set(gca,'ylim',[-0.04 0.04]);
legend('estimation','true');
title('time (s)','fontsize',13)

% %% with LPFilter
% figure;
% subplot(311);hold on;plot(Time,faults_in(1,1:k),'r--','linewidth',1.5);plot(Time,f_ob_filter(1,1:k),'b','linewidth',1.5);
% ylabel('f_{Ax} (m/s^2)','fontsize',13);grid;
% title('Accelerometer Faults Estimation via NDO and LowPassFilter')
% set(gca,'xlim',[0 d_t*k],'fontsize',13); %set(gca,'ylim',[-2 4]);
% subplot(312);hold on;plot(Time,faults_in(2,1:k),'r--','linewidth',1.5);plot(Time,f_ob_filter(2,1:k),'b','linewidth',1.5); 
% ylabel('f_{Ay} (m/s^2)','fontsize',13);grid;
% set(gca,'xlim',[0 d_t*k],'fontsize',13); %set(gca,'ylim',[-2 4]);
% subplot(313);hold on;plot(Time,faults_in(3,1:k),'r--','linewidth',1.5);plot(Time,f_ob_filter(3,1:k),'b','linewidth',1.5); 
% ylabel('f_{Az} (m/s^2)','fontsize',13);grid;
% set(gca,'xlim',[0 d_t*k],'fontsize',13); %set(gca,'ylim',[5000 8000]);
% legend('true','estimation');
% h1=axes('position',[0.52 0.0001 0.0001 0.0001],'fontsize',13);axis off;
% title('time (s)','fontsize',13)   
%     
% figure;
% subplot(311);hold on;plot(Time,faults_in(4,1:k),'r--','linewidth',1.5);plot(Time,f_ob_filter(4,1:k),'b','linewidth',1.5);
% ylabel('f_{p} (rad/s)','fontsize',13);grid;
% title('Gyroscope Faults Estimation via NDO and LowPassFilter')
% set(gca,'xlim',[0 d_t*k],'fontsize',13); % set(gca,'ylim',[-0.02 0.04]);
% subplot(312);hold on;plot(Time,faults_in(5,1:k),'r--','linewidth',1.5);plot(Time,f_ob_filter(5,1:k),'b','linewidth',1.5);
% ylabel('f_{q} (rad/s)','fontsize',13);grid;
% set(gca,'xlim',[0 d_t*k],'fontsize',13); % set(gca,'ylim',[5000 8000]);
% subplot(313);hold on;plot(Time,faults_in(6,1:k),'r--','linewidth',1.5);plot(Time,f_ob_filter(6,1:k),'b','linewidth',1.5);
% ylabel('f_{r} (rad/s)','fontsize',13);grid;
% set(gca,'xlim',[0 d_t*k],'fontsize',13);  %set(gca,'ylim',[-0.02 0.06]);
% legend('true','estimation');
% h1=axes('position',[0.52 0.0001 0.0001 0.0001],'fontsize',13);axis off;
% title('time (s)','fontsize',13)

%% END