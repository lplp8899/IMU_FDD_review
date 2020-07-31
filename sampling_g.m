
function [u_n_m] = sampling_g(u_n_m0,dt,dt_g)


% u_n_m0 = 
% dt = 0.01;
% dt_g = 0.01;
gen = length(u_n_m0);
u_n_m = zeros(1,gen);

n_scale = dt_g/dt;

% [N,res] = quorem(sym(gen),sym(n_scale));
N = gen/n_scale;
res = gen - N*n_scale;

for i = 1:N
    u_n_m(:,(i-1)*n_scale+1 : (i-1)*n_scale+n_scale) = u_n_m0(1, (i-1)*n_scale+1);
end



u_n_m(gen-res:gen) = u_n_m(gen-res);



% close all
% k = gen;
% Time = dt*(1:k);
% figure;
% subplot(211); hold on; plot(Time,u_n_m(1,1:k),'r','linewidth',1.5); plot(Time,u_n_m0(1,1:k),'b--','linewidth',1.5); ylabel('u_n (m/s)','fontsize',13);grid;
%     set(gca,'xlim',[0 delta_t*k],'fontsize',13);
% legend('true','estimation');
% subplot(212); hold on; plot(Time,u_n_m(1,1:k) - u_n_m0(1,1:k),'b--','linewidth',1.5); ylabel('u_n (m/s)','fontsize',13);grid;
%     set(gca,'xlim',[0 delta_t*k],'fontsize',13);
% % h1=axes('position',[0.05 0.05 0.89 0.89],'fontsize',12);axis off;title(h1,'States')
% h1=axes('position',[0.35 0.0001 0.0001 0.0001],'fontsize',13);axis off;
% title('time (s)','fontsize',13)
% h1=axes('position',[0.79 0.0001 0.0001 0.0001],'fontsize',13);
% title('time (s)','fontsize',13)









