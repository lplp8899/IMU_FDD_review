function y_ob = model_output(x,u)
% do not use the lat,lon,h as states

[num,~] = size(x);
% the real measurements
% lat = x(1,:); lon = x(2,:); h = x(3,:);
UN = x(4-3,:); UE = x(5-3,:); UD = x(6-3,:);
phi = x(7-3,:); theta = x(8-3,:); psi = x(9-3,:);
if num  == 15-3 || num  == 12-3
    la_x = x(10-3,:); la_y = x(11-3,:); la_z = x(12-3,:);
    la_p = x(13-3,:); la_q = x(14-3,:); la_r = x(15-3,:);
%     uwind = x(13,:);vwind = x(14,:);wwind = x(15,:);
elseif num  == 9-3
    la_x = 0; la_y = 0; la_z = 0;
    la_p = 0; la_q = 0; la_r = 0;
% elseif num  == 6
%     la_x = 0; la_y = 0; la_z = 0;
%     la_p = 0; la_q = 0; la_r = 0;
end
 

%];
 
y_ob = [
     UN;
     UE;
     UD;
     phi;
     theta;
     psi;];

 
 
 
 
 
 
 
 
       