% subfunctions 
% augmented states of the observer
function x_ob=model_fx(x,u,g,delta_t,lat,lon,h)

[num,~]  =  size(x);

UN = x(4-3,:); UE = x(5-3,:); UD = x(6-3,:);
phi = x(7-3,:); theta = x(8-3,:); psi = x(9-3,:);

if num  == 15 || num  == 12
    la_x = x(10-3,:); la_y = x(11-3,:); la_z = x(12-3,:);
    la_p = x(13-3,:); la_q = x(14-3,:); la_r = x(15-3,:);
%     uwind = x(13,:);vwind = x(14,:);wwind = x(15,:);
elseif num  == 9
    la_x = 0; la_y = 0; la_z = 0;
    la_p = 0; la_q = 0; la_r = 0;
elseif num  == 6
    la_x = 0; la_y = 0; la_z = 0;
    la_p = 0; la_q = 0; la_r = 0;
end
% input
Axm = u(1,:); Aym = u(2,:); Azm = u(3,:);
pm = u(4,:); qm = u(5,:); rm = u(6,:); 

% para
Omega = 7.2921e-5;
RN = 6375000;
RE = 6391000;
Re = 6367434;


%
if num==15
x_ob=[   ( (Axm - la_x).*cos(theta).*cos(psi) + (Aym - la_y).* ( sin(phi).*sin(theta).*cos(psi) - cos(phi).*sin(psi) ) ...
            + ( Azm - la_z ).*( cos(phi).*sin(theta).*cos(psi) + sin(phi).*sin(psi) ) );
         ( (Axm - la_x).*cos(theta).*sin(psi) + (Aym - la_y).* ( sin(phi).*sin(theta).*sin(psi) + cos(phi).*cos(psi) ) ...
            + ( Azm - la_z ).*( cos(phi).*sin(theta).*sin(psi) - sin(phi).*cos(psi) )  );
         ( - (Axm - la_x).*sin(theta) + (Aym- la_y).*sin(phi).*cos(theta) + (Azm- la_z).*cos(phi).*cos(theta)  + g );
         ( (pm - la_p) + (qm - la_q).*sin(phi).*tan(theta) + (rm - la_r).*cos(phi).*tan(theta) ...
            - ( UE./(Re + h) + Omega.*cos(lat) ).*cos(psi)./cos(theta) + UN.*sin(psi)./ ((Re + h).* cos(theta) ) );
         ( (qm - la_q).*cos(phi) - (rm - la_r).*sin(phi) + ( UE./(Re + h) + Omega.*cos(lat) ).*sin(psi) + UN.*cos(psi)./(Re + h) );% not r*sin(theta)
         ( (qm - la_q).*sin(phi)./cos(theta) + (rm - la_r).*cos(phi)./cos(theta)...
            + ( UE./(Re + h) + Omega*cos(lat) ).*tan(theta).*cos(psi) + UN.*tan(theta).*sin(psi)./(Re + h) + UE.*tan(lat)./(Re + h) + Omega.*sin(lat) );
%          ( (pm - la_p) + (qm - la_q).*sin(phi).*tan(theta) + (rm - la_r).*cos(phi).*tan(theta) );
%          ( (qm - la_q).*cos(phi) - (rm - la_r).*sin(phi)) ; % not r*sin(theta)
%          ( (qm - la_q).*sin(phi)./cos(theta) + (rm - la_r).*cos(phi)./cos(theta)) ;	
         0*ones(1,den); 
         0*ones(1,den); 
         0*ones(1,den); 
         0*ones(1,den); 
         0*ones(1,den); 
         0*ones(1,den);
         0*ones(1,den); 
         0*ones(1,den); 
         0*ones(1,den);];
elseif num==12
    x_ob=[ ( (Axm - la_x).*cos(theta).*cos(psi) + (Aym - la_y).* ( sin(phi).*sin(theta).*cos(psi) - cos(phi).*sin(psi) ) ...
            + ( Azm - la_z ).*( cos(phi).*sin(theta).*cos(psi) + sin(phi).*sin(psi) ) );
         ( (Axm - la_x).*cos(theta).*sin(psi) + (Aym - la_y).* ( sin(phi).*sin(theta).*sin(psi) + cos(phi).*cos(psi) ) ...
            + ( Azm - la_z ).*( cos(phi).*sin(theta).*sin(psi) - sin(phi).*cos(psi) )  );
         ( - (Axm - la_x).*sin(theta) + (Aym- la_y).*sin(phi).*cos(theta) + (Azm- la_z).*cos(phi).*cos(theta)  + g );
         ( (pm - la_p) + (qm - la_q).*sin(phi).*tan(theta) + (rm - la_r).*cos(phi).*tan(theta) ...
            - ( UE./(Re + h) + Omega.*cos(lat) ).*cos(psi)./cos(theta) + UN.*sin(psi)./ ((Re + h).* cos(theta) ) );
         ( (qm - la_q).*cos(phi) - (rm - la_r).*sin(phi) + ( UE./(Re + h) + Omega.*cos(lat) ).*sin(psi) + UN.*cos(psi)./(Re + h) );% not r*sin(theta)
         ( (qm - la_q).*sin(phi)./cos(theta) + (rm - la_r).*cos(phi)./cos(theta)...
            + ( UE./(Re + h) + Omega*cos(lat) ).*tan(theta).*cos(psi) + UN.*tan(theta).*sin(psi)./(Re + h) + UE.*tan(lat)./(Re + h) + Omega.*sin(lat) );
%          ( (pm - la_p) + (qm - la_q).*sin(phi).*tan(theta) + (rm - la_r).*cos(phi).*tan(theta) );
%          ( (qm - la_q).*cos(phi) - (rm - la_r).*sin(phi)) ; % not r*sin(theta)
%          ( (qm - la_q).*sin(phi)./cos(theta) + (rm - la_r).*cos(phi)./cos(theta)) ;	
         0*ones(1,den); 
         0*ones(1,den); 
         0*ones(1,den); 
         0*ones(1,den); 
         0*ones(1,den); 
         0*ones(1,den);];
elseif num==9
    x_ob=[   ( (Axm - la_x).*cos(theta).*cos(psi) + (Aym - la_y).* ( sin(phi).*sin(theta).*cos(psi) - cos(phi).*sin(psi) ) ...
            + ( Azm - la_z ).*( cos(phi).*sin(theta).*cos(psi) + sin(phi).*sin(psi) ) );
         ( (Axm - la_x).*cos(theta).*sin(psi) + (Aym - la_y).* ( sin(phi).*sin(theta).*sin(psi) + cos(phi).*cos(psi) ) ...
            + ( Azm - la_z ).*( cos(phi).*sin(theta).*sin(psi) - sin(phi).*cos(psi) )  );
         ( - (Axm - la_x).*sin(theta) + (Aym- la_y).*sin(phi).*cos(theta) + (Azm- la_z).*cos(phi).*cos(theta)  + g );
         ( (pm - la_p) + (qm - la_q).*sin(phi).*tan(theta) + (rm - la_r).*cos(phi).*tan(theta) ...
            - ( UE./(Re + h) + Omega.*cos(lat) ).*cos(psi)./cos(theta) + UN.*sin(psi)./ ((Re + h).* cos(theta) ) );
         ( (qm - la_q).*cos(phi) - (rm - la_r).*sin(phi) + ( UE./(Re + h) + Omega.*cos(lat) ).*sin(psi) + UN.*cos(psi)./(Re + h) );% not r*sin(theta)
         ( (qm - la_q).*sin(phi)./cos(theta) + (rm - la_r).*cos(phi)./cos(theta)...
            + ( UE./(Re + h) + Omega*cos(lat) ).*tan(theta).*cos(psi) + UN.*tan(theta).*sin(psi)./(Re + h) + UE.*tan(lat)./(Re + h) + Omega.*sin(lat) );
%          ( (pm - la_p) + (qm - la_q).*sin(phi).*tan(theta) + (rm - la_r).*cos(phi).*tan(theta) );
%          ( (qm - la_q).*cos(phi) - (rm - la_r).*sin(phi)) ; % not r*sin(theta)
%          ( (qm - la_q).*sin(phi)./cos(theta) + (rm - la_r).*cos(phi)./cos(theta)) ;	
         0*la_V; 
         0*la_a; 
         0*la_b;];
elseif num==6
    % update 1s
%     n = (k-1)/100;
%     nn = (k-1)/100;
%     if abs(n-round(n))<1e-6
%         x_ob=[   (Axm-la_x)-g.*sin(theta)+(rm-la_r).*vGS_b-(qm-la_q).*wGS_b;
%              (Aym-la_y)+g.*cos(theta).*sin(phi)+(pm-la_p).*wGS_b-(rm-la_r).*uGS_b;
%              (Azm-la_z)+g.*cos(theta).*cos(phi)+(qm-la_q).*uGS_b-(pm-la_p).*vGS_b;
%              (pm-la_p)+(qm-la_q).*sin(phi).*tan(theta)+(rm-la_r).*cos(phi).*tan(theta);
%              (qm-la_q).*cos(phi)-(rm-la_r).*sin(phi);
%              (qm-la_q).*sin(phi)./cos(theta)+(rm-la_r).*cos(phi)./cos(theta);];
%     elseif abs(nn-round(nn))<1e-6
%         x_ob=[   uGS_b;
%              (Aym-la_y)+g.*cos(theta).*sin(phi)+(pm-la_p).*wGS_b-(rm-la_r).*uGS_b;
%              wGS_b;
%              (pm-la_p)+(qm-la_q).*sin(phi).*tan(theta)+(rm-la_r).*cos(phi).*tan(theta);
%              (qm-la_q).*cos(phi)-(rm-la_r).*sin(phi);
%              (qm-la_q).*sin(phi)./cos(theta)+(rm-la_r).*cos(phi)./cos(theta);];
%     else
%         x_ob=[   uGS_b;
%              vGS_b;
%              wGS_b;
%              (pm-la_p)+(qm-la_q).*sin(phi).*tan(theta)+(rm-la_r).*cos(phi).*tan(theta);
%              (qm-la_q).*cos(phi)-(rm-la_r).*sin(phi);
%              (qm-la_q).*sin(phi)./cos(theta)+(rm-la_r).*cos(phi)./cos(theta);];
%     end
    %
    x_ob=[ ( (Axm - la_x).*cos(theta).*cos(psi) + (Aym - la_y).* ( sin(phi).*sin(theta).*cos(psi) - cos(phi).*sin(psi) ) ...
            + ( Azm - la_z ).*( cos(phi).*sin(theta).*cos(psi) + sin(phi).*sin(psi) ) );
         ( (Axm - la_x).*cos(theta).*sin(psi) + (Aym - la_y).* ( sin(phi).*sin(theta).*sin(psi) + cos(phi).*cos(psi) ) ...
            + ( Azm - la_z ).*( cos(phi).*sin(theta).*sin(psi) - sin(phi).*cos(psi) )  );
         ( - (Axm - la_x).*sin(theta) + (Aym- la_y).*sin(phi).*cos(theta) + (Azm- la_z).*cos(phi).*cos(theta)  + g );
         ( (pm - la_p) + (qm - la_q).*sin(phi).*tan(theta) + (rm - la_r).*cos(phi).*tan(theta) ...
            - ( UE./(Re + h) + Omega.*cos(lat) ).*cos(psi)./cos(theta) + UN.*sin(psi)./ ((Re + h).* cos(theta) ) );
         ( (qm - la_q).*cos(phi) - (rm - la_r).*sin(phi) + ( UE./(Re + h) + Omega.*cos(lat) ).*sin(psi) + UN.*cos(psi)./(Re + h) );% not r*sin(theta)
         ( (qm - la_q).*sin(phi)./cos(theta) + (rm - la_r).*cos(phi)./cos(theta)...
            + ( UE./(Re + h) + Omega*cos(lat) ).*tan(theta).*cos(psi) + UN.*tan(theta).*sin(psi)./(Re + h) + UE.*tan(lat)./(Re + h) + Omega.*sin(lat) ); ];
%          ( (pm - la_p) + (qm - la_q).*sin(phi).*tan(theta) + (rm - la_r).*cos(phi).*tan(theta) );
%          ( (qm - la_q).*cos(phi) - (rm - la_r).*sin(phi)) ; % not r*sin(theta)
%          ( (qm - la_q).*sin(phi)./cos(theta) + (rm - la_r).*cos(phi)./cos(theta)) ;	;];	
         
end