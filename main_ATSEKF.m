%-----------------------------------------------------------------------
%   This is the code for the following paper:
%
%   Q. He, W. Zhang, P. Lu and J. Liu
%   Performance comparison of representative model-based fault 
%   reconstruction algorithms for aircraft sensor fault detection 
%   and diagnosis
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

r2d = 180/pi;
delta_t = 0.01;
GPS_rate = [1; 10; 100];
num = 3;
S_hz = GPS_rate(num);
load new_sim_Ma_06_H_1000_notur.mat

%% organize data

% true states
Vtas   =  Out_nonoise(:,1)';
alpha  =  Out_nonoise(:,2)';
beta   =  Out_nonoise(:,3)';
pb     =  Out_nonoise(:,4)';
qb     =  Out_nonoise(:,5)';
rb     =  Out_nonoise(:,6)';
phi    =  Out_nonoise(:,9)';
theta  =  Out_nonoise(:,8)';
psi    =  Out_nonoise(:,7)';
xe     =  Out_nonoise(:,10)';
ye     =  Out_nonoise(:,11)';
ze     =  Out_nonoise(:,12)';
u_b    =  Out_nonoise(:,13)';
v_b    =  Out_nonoise(:,14)';
w_b    =  Out_nonoise(:,15)';
u_n    =  Out_nonoise(:,16)';
v_n    =  Out_nonoise(:,17)';
w_n    =  Out_nonoise(:,18)';

% Acc measurement with noise
A_x    =  squeeze(Ax_m)';
A_y    =  squeeze(Ay_m)';
A_z    =  squeeze(Az_m)';

% turbulence
ug_b     = Tur_b(:,1)';
vg_b     = Tur_b(:,2)';
wg_b     = Tur_b(:,3)';
% p_tur    = Tur_b(:,4)';

% derivative of the turbulence
ug_b_dot     = Tur_b_dot(:,1)';
vg_b_dot     = Tur_b_dot(:,2)';
wg_b_dot     = Tur_b_dot(:,3)';

% NOTE the derivative comes from last step!!! 
% otherwise will cause a delay
ug_b_dot      =  [0 ug_b_dot(1:length(ug_b_dot)-1)];
vg_b_dot      =  [0 vg_b_dot(1:length(ug_b_dot)-1)];
wg_b_dot      =  [0 wg_b_dot(1:length(ug_b_dot)-1)];

% ground velocity in the body axis 
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
p_m         =  Out_noise(:,4)';
q_m         =  Out_noise(:,5)';
r_m         =  Out_noise(:,6)';
Ax_m        =  squeeze(Ax_m)';
Ay_m        =  squeeze(Ay_m)';
Az_m        =  squeeze(Az_m)';
Vt_m        =  Out_noise(:,1)';
alpha_m     =  Out_noise(:,2)';
beta_m      =  Out_noise(:,3)';
phi_m       =  Out_noise(:,9)';
theta_m     =  Out_noise(:,8)';
psi_m       =  Out_noise(:,7)';
ze_m        =  Out_noise(:,12)';
xe_m        =  Out_noise(:,10)';
ye_m        =  Out_noise(:,11)';
% u_b_m       =  Out_noise(:,13)'; % not realistic
% v_b_m       =  Out_noise(:,14)'; % not realistic
% w_b_m       =  Out_noise(:,15)'; % not realistic
u_n_m       =  Out_noise(:,16)';
v_n_m       =  Out_noise(:,17)';
w_n_m       =  Out_noise(:,18)';

% different sampling rate
u_n_m0 = u_n_m;
v_n_m0 = v_n_m;
w_n_m0 = w_n_m;
dtu_g = 1/S_hz;
dtv_g = 1/S_hz;
dtw_g = 1/S_hz;
u_n_m = sampling_g(u_n_m0,delta_t,dtu_g);
v_n_m = sampling_g(v_n_m0,delta_t,dtv_g);
w_n_m = sampling_g(w_n_m0,delta_t,dtw_g);
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
x_real = [ u_n; v_n; w_n; phi; theta; psi;];

u_origin = u;
faults_in = zeros(6,gen);
fy = 0.001;

faults_in(1,1001:3000) = faults_in(1,1001:3000) + 1 * 2;
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

%% initialization of ATSEKF

z_real  =  [ UN_m; UE_m; UD_m; phi_m; theta_m; psi_m;];
x_ob_0 = [ UN_m(1); UE_m(1); UD_m(1); phi_m(1); theta_m(1); psi_m(1); ];

P_x0 = 1*eye(size(x_ob_0,1));
[dim_sys,~] = size(x_ob_0);
[dim_out,~] = size(z_real);
x_ob = zeros(dim_sys,2*dim_sys+1);
z_ob = zeros(dim_out,2*dim_sys+1);
x_ob_mean = zeros(dim_sys,gen);
x_ob_mean1 = zeros(dim_sys,gen);
z_ob_mean = zeros(dim_out,gen);
inno = zeros(dim_out,gen);
norm_inno = zeros(1,gen);
x_ob_filter = zeros(dim_sys,gen);
z_ob_filter = zeros(dim_out,gen);
P_ob_filter = zeros(dim_sys,dim_sys,gen);
residual = zeros(dim_out,gen);
norm_residual = zeros(1,gen);
error_x = zeros(6,gen);
norm_error_x = zeros(1,gen);
error_b = zeros(dim_u,gen);
S_inno = zeros(dim_out,dim_out,gen);
C_inno = zeros(dim_out,dim_out,gen);
norm_S_inno = zeros(1,gen); 
norm_C_inno = zeros(1,gen); 
S_resi = zeros(dim_out,dim_out,gen);
C_resi = zeros(dim_out,dim_out,gen);
norm_S_resi = zeros(1,gen); 
norm_C_resi = zeros(1,gen);
Qd_est = zeros(dim_sys,dim_sys,gen);
Qk_est = zeros(dim_sys,dim_sys,gen);
Qk = diag([1e-4,1e-4,1e-4,3e-8,3e-8,3e-8]);
Rk = diag([1e-4,1e-4,1e-4,3e-8,3e-8,3e-8]); 
max_iter=500;% record the iteration no. of IEKF
epsu_crit=1e-8;

epsu=zeros(1,5000);
adaQ1=zeros(1,5000);
adaQ2=zeros(1,5000);
adaQ3=zeros(1,5000);
adaQ4=zeros(1,5000);
adaQ5=zeros(1,5000);
adaQ6=zeros(1,5000);

dim_bias = 6;
b_ob=zeros(dim_u,gen);
Pb_ob=zeros(dim_u,dim_u,gen);
x_hat = zeros(dim_sys,1);
Px_hat = zeros(dim_sys,dim_sys,gen);
b_ob_0 = [ 1e-3; 1e-3; 1e-3; 1e-3; 1e-3; 1e-3]; % input faults
Pb0 = 1e0*eye(dim_u);
Pxb0 = 1e-6*eye(dim_sys);
V0 = Pxb0/Pb0;
x_ob_0 = x_ob_0-V0*b_ob_0;
Px_0 = P_x0-V0*Pb0*V0';
Pb_0 = Pb0;

% initialization of Q_b and Q_xb
%Q_b = diag([0.01^1,0.01^2,0.01^3,(0.01/57.3)^3,(0.01/57.3)^1,(0.01/57.3)^2]);
%Q_b = diag([0.01^2,0.01^2,0.01^2,(0.01/57.3)^2,(0.01/57.3)^2,(0.01/57.3)^2]);
%Q_b = diag([0.02^3,0.01^3,0.01^3,(0.02/57.3)^3,(0.02/57.3)^3,(0.02/57.3)^3]);
%Q_b = diag([0.01^1,0.01^2,0.01^3,(0.01/57.3)^3,(0.01/57.3)^1,(0.01/57.3)^2]);
Q_b = diag([0.01^2,0.01^2,0.01^2,(0.01/57.3)^2,(0.01/57.3)^2,(0.01/57.3)^2]);
Q_xb = 0;

lat_g = 0;
lon_g = 0;
h_g = 0;
g = 9.80665;

%% ATSEKF

for k=1:5000
    % continuous state transition matrix
    F=F_Ja(x_ob_0,u(:,k),g);
    
    % continuous noise distribution matrix
    G_noise=G_Ja(x_ob_0); 
    
    % discrete matrices
    [Phi,Gamma]=c2d(F,G_noise,delta_t);
 
    % bias term in the bias-free filter
    F_in = Gamma;
    U_bar = Phi*V0+F_in;
    Pb_k0 = Pb_0 + Q_b;
    U = U_bar + (Q_xb-U_bar*Q_b)/(Pb_k0);
    u_0 = (U_bar-U)*b_ob_0;
    Q_bar = Gamma*Qk*Gamma' - Q_xb*U_bar'-U*(Q_xb-U_bar*Q_b)';   
    P =Phi*P_x0*Phi'+ Q_bar;
    
    % predict
    x_ob = x_ob_0 + model_fx(x_ob_0,u(:,k),g,d_t,lat_g,lon_g,h_g)*delta_t + u_0;
    
    % iteration
    eta_1=x_ob;
    num_iter=0;
    flag=1;
    while flag==1
        H=H_Ja(eta_1);        
        V=H*P*H'+Rk;
        K=P*H'/V;
        z_ob=model_output(eta_1,u);

        % residual
        inno(:,k)=z_real(:,k)-z_ob;
        norm_inno(:,k)=norm(inno(:,k));

        eta_2=x_ob+K*(inno(:,k)-H*(x_ob-eta_1));
        epsu(:,k)=norm(eta_2-eta_1)./norm(eta_2);

        if (k<=20)&&(epsu(:,k)>epsu_crit) && (num_iter<max_iter)
            eta_1=eta_2;
            num_iter=num_iter+1;
            disp('heehee')
        else 
            flag=0;
        end
    end
   
    x_ob_iekf=eta_2;
    P_ob_iekf=(eye(dim_sys)-K*H)*P*(eye(dim_sys)-K*H)'+K*Rk*K';
    L = K;
    x_ob_filter(:,k) = x_ob_iekf;
    P_ob_filter(:,:,k) = P_ob_iekf;
        
    % bias filter
    S = H*U;
    Kb = Pb_k0*S'/(H*P*H'+Rk+S*Pb_k0*S');
    b_ob(:,k) = b_ob_0+Kb*(inno(:,k)-S*b_ob_0);
    Pb_ob(:,:,k) = (eye(dim_bias)-Kb*S)*Pb_k0;
    
    % estimation error of faults
    error_b(:,k) = faults_in(:,k) - b_ob(:,k);
  
    % coupled equations    
    V = U - L*S;
    x_hat(:,k) = x_ob_filter(:,k)+V*b_ob(:,k);
    Px_hat(:,:,k) = P_ob_filter(:,:,k)+V*Pb_ob(:,:,k)*V';
    
    % next generation
    x_ob_0=x_ob_filter(:,k);
    P_x0=P_ob_filter(:,:,k);
    b_ob_0 = b_ob(:,k);
    Pb_0 = Pb_ob(:,:,k);
    V0 = V;
    
    % new residual
    z_ob_filter(:,k)=model_output(x_hat(:,k),u);
    residual(:,k)=z_real(:,k)-z_ob_filter(:,k);
    norm_residual(:,k)=norm(residual(:,k));
    
    % Q adapatation 
    if k > 1
        inno(:,k) = z_real(:,k) - x_hat(:,k-1) - model_fx( (x_hat(:,k-1) ),u(:,k),g,d_t,lat_g,lon_g,h_g)*delta_t - Gamma*b_ob(:,k-1)*delta_t;
        inno(:,k) = z_real(:,k) - x_hat(:,k-1) - model_fx( (x_hat(:,k-1) ),u(:,k),g,d_t,lat_g,lon_g,h_g)*delta_t - Gamma*b_ob(:,k-1);
        inno(:,k) = z_real(:,k) - Phi*x_hat(:,k-1) - Gamma*b_ob(:,k-1); % approximated
    else
        inno(:,k) = 0;
    end

    S_inno(:,:,k) = inno(:,k)*inno(:,k)';
    S_resi(:,:,k) = residual(:,k)*residual(:,k)';
    N_win = 10;

    if N_win == 1
        C_inno(:,:,k) = S_inno(:,:,k);
        C_resi(:,:,k) = S_resi(:,:,k);
    end

    if N_win > 1
        if k == 1
            C_inno(:,:,k) = S_inno(:,:,k);
            C_resi(:,:,k) = S_resi(:,:,k);
        elseif k >1 && k < N_win
            C_inno(:,:,k) = sum(S_inno(:,:,1:k),3)/(k-1);
            C_resi(:,:,k) = sum(S_resi(:,:,1:k),3)/(k-1);
        elseif k >= N_win
            C_inno(:,:,k) = sum(S_inno(:,:,k-N_win+1:k),3)/(N_win-1);
            C_resi(:,:,k) = sum(S_resi(:,:,k-N_win+1:k),3)/(N_win-1);
        end 
    end

    norm_S_inno(:,k) = norm(S_inno(:,:,k));
    norm_C_inno(:,k) = norm(C_inno(:,:,k));
    norm_S_resi(:,k) = norm(S_resi(:,:,k));
    norm_C_resi(:,k) = norm(C_resi(:,:,k));
    num_G = 0; 

    if num_G == 1
        Gamma_H = H * G_noise * 1; 
        Qd_0 = C_resi(:,:,k) - Gamma_H*Qk*Gamma_H'*delta_t - Rk;

        Qd_est = diag([ max(Qd_0(1,1),0), max(Qd_0(2,2),0), max(Qd_0(3,3),0),...
                    max(Qd_0(4,4),0), max(Qd_0(5,5),0), max(Qd_0(6,6),0)]);

        inv_Gamma = eye(dim_u)/Gamma_H;
        Qdx_est = inv_Gamma * Qd_est *inv_Gamma' / (delta_t^1);
        Qk_est(:,:,k+1) = diag([ Qdx_est(1,1), Qdx_est(2,2), Qdx_est(3,3), ...
                                Qdx_est(4,4), Qdx_est(5,5), Qdx_est(6,6)]); 

    else
        Gamma_H = H * G_noise * delta_t; 
        Qd_0 = C_resi(:,:,k) - Gamma_H*Qk*Gamma_H'/delta_t - Rk;
        Qd_est = diag([ max(Qd_0(1,1),0), max(Qd_0(2,2),0), max(Qd_0(3,3),0),...
                       max(Qd_0(4,4),0), max(Qd_0(5,5),0), max(Qd_0(6,6),0)]);
        inv_Gamma = eye(dim_u)/Gamma_H;
        Qdx_est = inv_Gamma * Qd_est *inv_Gamma' * (delta_t^2);
        Qk_est(:,:,k+1) = diag([ Qdx_est(1,1), Qdx_est(2,2), Qdx_est(3,3), ...
                             Qdx_est(4,4), Qdx_est(5,5), Qdx_est(6,6)]); 
    end

    Q_b = Qk_est(:,:,k+1);

end
toc;

%% save the adaptive Q_b for ploting

for i=1:5000
    adaQ1(i)=Qk_est(1,1,i);
    adaQ2(i)=Qk_est(2,2,i);
    adaQ3(i)=Qk_est(3,3,i);
    adaQ4(i)=Qk_est(4,4,i);
    adaQ5(i)=Qk_est(5,5,i);
    adaQ6(i)=Qk_est(6,6,i);    
end

%% plot 

Time=delta_t*(1:k);

figure;
subplot(311);hold on;plot(Time,b_ob(1,1:k),'b','linewidth',1.5);plot(Time,faults_in(1,1:k),'r--','linewidth',1.5); 
ylabel('f_{Ax} (m/s^2)','fontsize',13);grid;
set(gca,'xlim',[0 delta_t*k],'fontsize',13);set(gca,'ylim',[-1.5 3.5]);
subplot(312);hold on;plot(Time,b_ob(2,1:k),'b','linewidth',1.5);plot(Time,faults_in(2,1:k),'r--','linewidth',1.5);
ylabel('f_{Ay} (m/s^2)','fontsize',13);grid;
set(gca,'xlim',[0 delta_t*k],'fontsize',13);set(gca,'ylim',[-1.5 2.5]);
subplot(313);hold on;plot(Time,b_ob(3,1:k),'b','linewidth',1.5);plot(Time,faults_in(3,1:k),'r--','linewidth',1.5);
ylabel('f_{Az} (m/s^2)','fontsize',13);grid;
set(gca,'xlim',[0 delta_t*k],'fontsize',13);set(gca,'ylim',[-2.2 2.2]);
legend('estimation','true');
title('time (s)','fontsize',13)    

figure;
subplot(311);hold on;plot(Time,b_ob(4,1:k),'b','linewidth',1.5);plot(Time,faults_in(4,1:k),'r--','linewidth',1.5);
ylabel('f_{p} (rad/s)','fontsize',13);grid;
set(gca,'xlim',[0 delta_t*k],'fontsize',13);set(gca,'ylim',[-0.03 0.06]);
subplot(312);hold on;plot(Time,b_ob(5,1:k),'b','linewidth',1.5);plot(Time,faults_in(5,1:k),'r--','linewidth',1.5);
ylabel('f_{q} (rad/s)','fontsize',13);grid;
set(gca,'xlim',[0 delta_t*k],'fontsize',13);set(gca,'ylim',[-0.03 0.045]);
subplot(313);hold on;plot(Time,b_ob(6,1:k),'b','linewidth',1.5);plot(Time,faults_in(6,1:k),'r--','linewidth',1.5);
ylabel('f_{r} (rad/s)','fontsize',13);grid;
set(gca,'xlim',[0 delta_t*k],'fontsize',13);set(gca,'ylim',[-0.04 0.04]);
legend('estimation','true');
title('time (s)','fontsize',13)  

figure;
subplot(321);plot(Time,adaQ1,'b','LineWidth',1.5);ylabel('Q_f of A_x channel');grid;
subplot(323);plot(Time,adaQ2,'b','LineWidth',1.5);ylabel('Q_f of A_y channel');grid;
subplot(325);plot(Time,adaQ3,'b','LineWidth',1.5);ylabel('Q_f of A_z channel');grid;
subplot(322);plot(Time,adaQ4,'b','LineWidth',1.5);ylabel('Q_f of p channel');grid;
subplot(324);plot(Time,adaQ5,'b','LineWidth',1.5);ylabel('Q_f of q channel');grid;
subplot(326);plot(Time,adaQ6,'b','LineWidth',1.5);ylabel('Q_f of r channel');grid;
%% END