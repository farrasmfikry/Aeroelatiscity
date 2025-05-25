%% Changing variables
clear

V_0 = 8;                     % Mean wind speed, constant for the axial induction factor
T = 300;                     % Total time [sec]
delta_t = 0.1;               % Time step [sec]

v = 0;                       % Wind Shear
turbulence_model = true;     % True or false (if false, then wind shear only)
tower_model = false;         % True or false (if false, no tower model)
dynamic_stall = true;        % True or false (if false, no dynamic stall)
pitch_controller = false;    % True or false (if false, no pitch controller and constant omega)
structural_response = true;  % True or false (if false, no structural response)
gravity = false;             % True or false (if false, no gravity)

theta_tilt = deg2rad(0);     % Tilt angle [rad]
theta_cone = deg2rad(0);     % Cone angle [rad]
theta_yaw  = deg2rad(0);     % Yaw angle [rad]
% theta_p = 0;               % Pitch angle [deg]
theta_initial = 0;

R    = 89.17;                % Blade radius [m]
H_WT = 119;                  % Hub height [m]
L_s  = 7.1;                  % Shaft length [m]
rho  = 1.225;

% pitch regulated controller (region 3)
K1 = 14;                  % Aerodynamic gain sechduling [deg]
K_P = 1.5;                % proportional gain in region 3 [rad/(rad/s)]
K_I = 0.64;               % intergral gain in region 3 [rad/rad]
omega_ref     = 1.05;     % controller reference rotational speed [rad/s] (need to update)
theta_range   = [0 45];   % physical stops for the pitch controller [deg]
dtheta_p_max  = 10;       % maximum pitch velocity [deg/s]
I = 1.6*10^8;             % inertia moment of the drivetrain [kg*m^2]

% optimal Cp tracking controller (region 1)
TSR = 8;                                     % design tsr = 8; design pitch = 0;
C_P = (10.64*10^6)/(0.5*rho*R^2*pi*11.4^3);  % according to 11.4 m/s as rated wind speed from handbook

K_opt = 0.5 * rho * R^2 * pi * R^3 * C_P / (TSR^3);
omega_rated = TSR*11.4/R;

omega_initial = 0.7;

M = 446000;     % mass of the nacelle, gearbox, generator [kg]      
k = 1.7*10^6;   % spring stiffness [N/m]

% alpha_filt = exp(-delta_t*pi/2);

%% Turbulence generator

% Number of points in box
n1=4096; % z-direction % 8192 in the case for 18 m/s and 4096 in the case for 7 m/s setup
n2=32;   % y-direction
n3=32;   % x-direction

% Size of box (make sure to change to your setup!)
if V_0 == 8
    Lz = 3276;
elseif V_0 == 18
    Lz = 7371;
else
    print("Wrong wind speed!")
end
Ly=180;
Lx=180;

% Change filename path if the folder structure is different
filename = "Turbulence_generator/sim1_v" + num2str(V_0) + ".bin";

%% Input Rotor Geometry

% Assignment 1 Data
bladedat = load("blade_data/bladedat.txt");
r = bladedat(:,1);      % Blade radius [m]
beta = bladedat(:,2);   % Beta twist angle [deg]
c = bladedat(:,3);      % Chord [m]
t = bladedat(:,4);      % Thickness-Chord-Ratio [%]
B = 3;                  % No. of blades [-]

r_add = r;

% Shorten (delete) the last element
r(end) = [];
c(end) = [];
t(end) = [];
beta(end) =[];

N = T/delta_t;              % Number of loops

theta_blade = zeros(3,N);   % Blade position (each column denotes a blade position) [rad]
theta_blade(1,1) = 0;       % Initial blade 1 position [rad]
theta_blade(2,1) = 2*pi/3;  % Initial blade 2 position [rad]
theta_blade(3,1) = 4*pi/3;  % Initial blade 3 position [rad]

% TRANSFORMATIONS
a_1 = [1 0 0;0 cos(theta_yaw) sin(theta_yaw);0 -sin(theta_yaw) cos(theta_yaw)];
a_2 = [cos(theta_tilt) 0 -sin(theta_tilt);0 1 0;sin(theta_tilt) 0 cos(theta_tilt)];

% Transformation from yaw+tilt (i.e. from ground to shaft)
a_12 = a_1*a_2;    

% Transformation from cone (i.e. from blade to cone)
a_34 = [cos(theta_cone) 0 -sin(theta_cone);0 1 0;sin(theta_cone) 0 cos(theta_cone)]; 

% The other way around
a_43 = transpose(a_34);
a_21 = transpose(a_12);

% Transformation from blade angle a_23 will be calculated in the time loop for
% each blade


%% Input Interpolation Constants

[aoa,cl,cd,cm,f_stat,cl_inv,cl_fs]=...
 readairfoildata_stat('blade_data/', ...
                    'cylinder.txt','FFA-W3-600.txt',...
                    'FFA-W3-480.txt','FFA-W3-360.txt',...
                    'FFA-W3-301.txt','FFA-W3-241.txt');
thick_prof = [100,60,48,36,30.1,24.1];

%% Load modeshapes data and structural response initialization 
data = readmatrix('blade_data/modeshapes.txt');

% GM1 = 1453
% GM2 = 2063
% GM3 = 660

u1fy = data(:,2);
u1fz = data(:,3);
u1ey = data(:,4);
u1ez = data(:,5);
u2fy = data(:,6);
u2fz = data(:,7);
m = data(:,8);

GM1 = trapz(r_add,u1fy.*u1fy.*m) + trapz(r_add,u1fz.*u1fz.*m);
GM2 = trapz(r_add,u1ey.*u1ey.*m) + trapz(r_add,u1ez.*u1ez.*m);
GM3 = trapz(r_add,u2fy.*u2fy.*m) + trapz(r_add,u2fz.*u2fz.*m);
GM = [GM1 GM2 GM3];

x_0 = zeros(5,1);
dx_0 = zeros(5,1);
dx_0(2,1) = omega_initial;

beta_newmark = 0.27;
gamma_newmark = 0.51;

tf = 1e-6;

massmatrix(1,1) = M+3*trapz(r_add,m);

massmatrix(3,1) = trapz(r_add,m.*u1fz);
massmatrix(4,1) = trapz(r_add,m.*u1ez);
massmatrix(5,1) = trapz(r_add,m.*u2fz);
massmatrix(1,3) = massmatrix(3,1);
massmatrix(1,4) = massmatrix(4,1);
massmatrix(1,5) = massmatrix(5,1);

massmatrix(2,3) = trapz(r_add,m.*r_add.*u1fy);
massmatrix(2,4) = trapz(r_add,m.*r_add.*u1ey);
massmatrix(2,5) = trapz(r_add,m.*r_add.*u2fy);
massmatrix(3,2) = massmatrix(2,3);
massmatrix(4,2) = massmatrix(2,4);
massmatrix(5,2) = massmatrix(2,5);

massmatrix(2,2) = trapz(r_add,r_add.*r_add.*m)*3;
massmatrix(3,3) = GM(1);
massmatrix(4,4) = GM(2);
massmatrix(5,5) = GM(3);

omega1f = 3.93;
omega1e = 6.10;
omega2f = 11.28;

stiffmatrix(1,1) = k;
stiffmatrix(2,2) = 0.01;
stiffmatrix(3,3) = omega1f^2*GM(1);
stiffmatrix(4,4) = omega1e^2*GM(2);
stiffmatrix(5,5) = omega2f^2*GM(3);

damp_factor = 0.03;
% dampmatrix(3,3) = omega1f*damp_factor*GM(1)/pi;
% dampmatrix(4,4) = omega1e*damp_factor*GM(2)/pi;
% dampmatrix(5,5) = omega2f*damp_factor*GM(3)/pi;
dampmatrix(1,1) = 0;
dampmatrix(2,2) = 0;
dampmatrix(3,3) = 0;
dampmatrix(4,4) = 0;
dampmatrix(5,5) = 0;

% Initial condition

ddx = zeros(N,5);
dx = zeros(N,5);
x = zeros(N,5);


%% Loop time

% Wind speed variation
Power_converged = zeros(length(V_0),1);
theta_p_converged = zeros(length(V_0),1);
omega_converged = zeros(length(V_0),1);
C_P_converged = zeros(length(V_0),1);

% time series
Power_Q2 = zeros(N,length(V_0));
omega_Q2 = zeros(N,length(V_0));
theta_p_Q2 = zeros(N,length(V_0));
Thrust_Q2  = zeros(N,length(V_0));
M_gen_Q2 = zeros(N,length(V_0));

TI = zeros(3,1);

for wsp = 1:length(V_0)
    
    time = zeros(N,1);
    r_1 = zeros(length(r),3*B,N);     % Position (x,y,z) in system 1 (column 1-3 blade 1, column 4-6 blade 2, column 7-9 blade 3)
    V_0_s4 = zeros(length(r),3*B,N);  % Velocity (V0_x,V0_y,V0_z) in system 4 (column 1-3 blade 1, column 4-6 blade 2, column 7-9 blade 3)
    p_z = zeros(length(r)+1,B,N);     % Normal Load
    p_y = zeros(length(r)+1,B,N);     % Tangential Load
    Power = zeros(N,1);               % Power
    Thrust = zeros(N,1);              % Thrust
    Thrust_tmp = zeros(N,B);          % Thrust on each blade
    theta_p = zeros(N+1,1);           % pitch angle
    omega = zeros(N+1,1);             % rotor speed
    omega(1) = omega_initial;
    x_tower = zeros(N+1,1);
    x_tower_dot = zeros(N+1,1);
    omega(1) = dx_0(2);         % Initial rotor speed
    theta_p_I = zeros(N+1,1);
    theta_p_P = zeros(N+1,1);
    M_gen = zeros(N+1,1);
    u_z_dot_dot = zeros(length(r_add),B,N);
    u_y_dot_dot = zeros(length(r_add),B,N);
    u_z_dot = zeros(length(r_add),B,N);
    u_y_dot = zeros(length(r_add),B,N);
    u_z = zeros(length(r_add),B,N);
    u_y = zeros(length(r_add),B,N);

    % theta_p(1) = theta_initial;

    % Induced wind initialization (W_y[1] and W_z[1] are initial values, representing n-1, therefore N+1 is needed)
    W_y_qs = zeros(length(r),B,N+1);
    W_z_qs = zeros(length(r),B,N+1);
    W_y_int = zeros(length(r),B,N+1);
    W_z_int = zeros(length(r),B,N+1);
    W_y = zeros(length(r),B,N+1);
    W_z = zeros(length(r),B,N+1);

    % Dynamic Stall
    f_s_matrix = zeros(3,length(r));
    alpha_store = zeros(length(r),N);

    % Turbulence box generator
    [X_turb,Y_turb,Z_turb,u] = turbulence_generator(filename,n1,n2,n3,Lx,Ly,Lz,H_WT,V_0);

    [X_turb_mesh,Y_turb_mesh] = meshgrid(X_turb,Y_turb); % Meshing the X,Y for 2D interpolation


    for n = 1:N
        time(n) = (n-1)*delta_t; % Loop through time

        if n > 1
            theta_blade(1,n) = theta_blade(1,n-1)+omega(n)*delta_t;
            theta_blade(2,n) = theta_blade(1,n-1)+2*pi/3;
            theta_blade(3,n) = theta_blade(1,n-1)+4*pi/3;
        end

        for b = 1:B

            % Transformations matrices for each blade
            a_23 = trans_a23(theta_blade(b,n));
            a_14 = a_34*a_23*a_12;
            a_41 = a_14';

            for i = 1:length(r)

                % (x,y,z) index for blades
                b_x_y_z = b*B-2:b*B;
                b_x = b*B-2;
                b_y = b*B-1;
                b_z = b*B;

                % Position in system 1
                r_1(i,b_x_y_z,n) = transpose([H_WT;0;0] + a_21*[0;0;-L_s] + a_41*[r(i);0;0]);

                % Turbulence model or only wind shear %
                if turbulence_model
                    % Turbulence model
                    u_plane = squeeze(u(n,:,:));
                    U_turb = interp2(X_turb_mesh,Y_turb_mesh,u_plane,r_1(i,b_x,n),r_1(i,b_y,n));
                    V_0_s1 = [0 0 V_0(wsp)*(r_1(i,b_x,n)/H_WT).^v+U_turb];
                else
                    % Wind Shear only
                    V_0_s1 = [0 0 V_0(wsp)*(r_1(i,b_x,n)/H_WT).^v];
                end

                if tower_model
                    % Wind Shear + Tower Influence for each blade
                    if r_1(i,b_x,n) > H_WT
                        a_r = 0;
                    else
                        a_r = 3.32;
                    end

                    r_t = sqrt(r_1(i,b_z,n)^2 + r_1(i,b_y,n)^2);

                    V_r =  (r_1(i,b_z,n)/r_t)*V_0_s1(3)*(1-(a_r/r_t)^2);
                    V_theta = (r_1(i,b_y,n)/r_t)*V_0_s1(3)*(1+(a_r/r_t)^2);

                    V_0_s1(2) = (r_1(i,b_y,n)/r_t)*V_r - (r_1(i,b_z,n)/r_t)*V_theta;
                    V_0_s1(3) = (r_1(i,b_z,n)/r_t)*V_r + (r_1(i,b_y,n)/r_t)*V_theta;

                end

                % Transition wind velocity from system 1 to system 4 for each blade
                V_0_s4(i,b_x_y_z,n) = transpose(a_14*V_0_s1');

                % BEM %
                V_rel_y = V_0_s4(i,b_y,n) + W_y(i,b,n) - omega(n)*r(i)*cos(theta_cone) - u_y_dot(i,b,n);
                V_rel_z = V_0_s4(i,b_z,n) + W_z(i,b,n) - u_z_dot(i,b,n) - x_tower_dot(n);

                phi = atand(V_rel_z./-V_rel_y);

                alpha = phi - (beta(i)+theta_p(n));

                V_rel_mag = sqrt((V_rel_y)^2+(V_rel_z)^2);

                if dynamic_stall
                    [~,C_d,~,fstat,clinv,clfs] = interpolation(aoa,cl,cd,cm,alpha,t(i),thick_prof,f_stat,cl_inv,cl_fs);

                    tau = 4*c(i)/V_rel_mag;

                    if n == 1  % Initialize
                        f_s_matrix(b,i) = fstat;
                    else       % Update
                        f_s_matrix(b,i) = fstat + (f_s_matrix(b,i)-fstat)*exp(-delta_t/tau);
                    end

                    C_l = f_s_matrix(b,i)*clinv + (1-f_s_matrix(b,i))*clfs;
                else
                    [C_l,C_d,~,~,~,~] = interpolation(aoa,cl,cd,cm,alpha,t(i),thick_prof,f_stat,cl_inv,cl_fs);
                end

                a = -W_z(i,b,n)/V_0(wsp) ; % Constant V_0

                if a <= 1/3
                    f_g = 1;
                else
                    f_g = (1/4)*(5-3*a);
                end

                V_0_W_n_mag =sqrt(V_0_s4(i,b_y,n).^2 +(V_0_s4(i,b_z,n)+f_g*W_z(i,b,n)).^2);

                l = 0.5*rho*V_rel_mag^2*c(i)*C_l;
                d = 0.5*rho*V_rel_mag^2*c(i)*C_d;

                p_z(i,b,n) = l*cos(deg2rad(phi))+d*sin(deg2rad(phi));
                p_y(i,b,n) = l*sin(deg2rad(phi))-d*cos(deg2rad(phi));

                W_z_qs(i,b,n+1) = -B*l*cos(deg2rad(phi))/(4*pi*rho*r(i)*Prandtl(B,R,r(i),deg2rad(phi))*V_0_W_n_mag);
                W_y_qs(i,b,n+1) = -B*l*sin(deg2rad(phi))/(4*pi*rho*r(i)*Prandtl(B,R,r(i),deg2rad(phi))*V_0_W_n_mag);

                % Dynamic wake/inflow model
                if a > 0.5
                    a_tau = 0.5;
                else
                    a_tau = a;
                end

                tau_1 = (1.1/(1-1.3*a_tau))*(R/V_0(wsp));
                tau_2 = (0.39-0.26*(r(i)/R)^2)*tau_1;
                k = 0.6;

                W_qs_current = [W_z_qs(i,b,n+1);W_y_qs(i,b,n+1)];
                W_qs_previous = [W_z_qs(i,b,n);W_y_qs(i,b,n)];

                H = W_qs_current + k.*tau_1.*((W_qs_current-W_qs_previous).*(1/delta_t));

                W_int_current = H + ([W_z_int(i,b,n);W_y_int(i,b,n)] - H).*exp(-delta_t/tau_1);
                W_z_int(i,b,n+1) = W_int_current(1);
                W_y_int(i,b,n+1) = W_int_current(2);

                W_previous = [W_z(i,b,n);W_y(i,b,n)];
                W_current = W_int_current + (W_previous-W_int_current).*exp(-delta_t/tau_2);
                W_z(i,b,n+1) = W_current(1);
                W_y(i,b,n+1) = W_current(2);

                alpha_store(i,n) = alpha;

                % BEM %

            end
        
            Power_tmp = omega(n)*trapz(r_add,p_y(:,b,n).*r_add);
            Power(n) = Power(n) + Power_tmp;
            Thrust_tmp(n,b) = trapz(r_add,p_z(:,b,n));
            Thrust(n) = Thrust(n) + Thrust_tmp(n,b);
            
          
        end

        disp([' n= ',num2str(n),', V= ',num2str(V_0(wsp))])

        if pitch_controller

            % region 3 pitch controller
            GK      = 1/(1 + theta_p(n)/K1);
            theta_p_P(n+1) = GK*K_P*(omega(n)-omega_ref);   % proportional error
            theta_p_I(n+1) = theta_p_I(n) + GK*K_I*(omega(n)- omega_ref)*delta_t;  % intergral error

            if theta_p_I(n+1) >= theta_range(2)
                theta_p_I(n+1) = theta_range(2);
            elseif theta_p_I(n+1) <= theta_range(1)
                theta_p_I(n+1) = theta_range(1);
            end

            theta_setpoint = theta_p_P(n+1) + theta_p_I(n+1);   % combine the error

            theta_p(n+1) = theta_setpoint;

            if theta_p(n+1) >= theta_p(n) + dtheta_p_max*delta_t
                theta_p(n+1) = theta_p(n) + dtheta_p_max*delta_t;
            elseif theta_p(n+1) < theta_p(n) - dtheta_p_max*delta_t
                theta_p(n+1) = theta_p(n) - dtheta_p_max*delta_t;
            end

            if theta_p(n+1) >= theta_range(2)
                theta_p(n+1) = theta_range(2);
            elseif theta_p(n+1) <= theta_range(1)
                theta_p(n+1) = theta_range(1);
            end
        end
        
        % Newmark method applied here
        if structural_response
          
            % add gravity force to the structural response

            if gravity
                p_y(:,1,n) = p_y(:,1,n) + m.*9.81.*sin(theta_blade(1,n));
                p_y(:,2,n) = p_y(:,2,n) + m.*9.81.*sin(theta_blade(2,n));
                p_y(:,3,n) = p_y(:,3,n) + m.*9.81.*sin(theta_blade(3,n));
            end

            if  omega(n) < omega_rated
                M_gen(n) = K_opt * omega(n)^2;
            else
                M_gen(n) = 10.64e6 / omega_rated;
            end

            if time(n) >= 100
                M_gen(n) = 0;
            end

            GF(1) = Thrust(n);
            GF(2) = (trapz(r_add,p_y(:,1,n).*r_add)+trapz(r_add,p_y(:,2,n).*r_add)+trapz(r_add,p_y(:,3,n).*r_add))-M_gen(n);
            GF(3) = trapz(r_add,p_y(:,1,n).*u1fy) + trapz(r_add,p_z(:,1,n).*u1fz);
            GF(4) = trapz(r_add,p_y(:,1,n).*u1ey) + trapz(r_add,p_z(:,1,n).*u1ez);
            GF(5) = trapz(r_add,p_y(:,1,n).*u2fy) + trapz(r_add,p_z(:,1,n).*u2fz);
            
            iter = 0;
            k = 1;

            if n == 1
                ddx_0 = massmatrix\(GF'-dampmatrix*dx_0-stiffmatrix*x_0);

                % Prediction step
                ddx(n,:) = ddx_0';
                dx(n,:) = dx_0'+delta_t.*ddx_0';
                x(n,:) = x_0'+delta_t.*dx_0'+0.5.*delta_t^2.*ddx_0';
            else
                % Prediction step
                ddx(n,:) = ddx(n-1,:);
                dx(n,:) = dx(n-1,:)+delta_t.*ddx(n-1,:);
                x(n,:) = x(n-1,:)+delta_t.*dx(n-1,:)+0.5.*delta_t^2.*ddx(n-1,:);
            end

            while k == 1
                iter = iter + 1;

                % Calculate residual
                residual = GF' - massmatrix*ddx(n,:)' - dampmatrix*dx(n,:)' - stiffmatrix*x(n,:)';
                r_max = max(abs(residual));

                % System matrices and increment correction
                K_star = stiffmatrix + ((gamma_newmark*delta_t)/(beta_newmark*delta_t^2)).*dampmatrix + (1/(beta_newmark*delta_t^2)).*massmatrix;

                partial_u = K_star\residual;

                x(n,:) = x(n,:)+partial_u';
                dx(n,:) = dx(n,:)+((gamma_newmark*delta_t)/(beta_newmark*delta_t^2))*partial_u';
                ddx(n,:) = ddx(n,:) + (1/(beta_newmark*delta_t^2))*partial_u';

                if r_max > tf
                    k = 1;
                else
                    k = 0;
                end

                if iter > 600
                    break
                end
            end

            x_tower(n+1) = x(n,1);
            x_tower_dot(n+1) = dx(n,1);
            omega(n+1) = dx(n,2);
            % omega(n+1) = omega(n) + (Power(n) - M_gen(n) - massmatrix(2,3)*ddx(n,3) ...
            %     - massmatrix(2,4)*ddx(n,4) - massmatrix(2,5)*ddx(n,5) ) / (I*M) * delta_t;

            u_y_dot_dot(:,1,n+1) = ddx(n,3).*u1fy + ddx(n,4).*u1ey + ddx(n,5).*u2fy;
            u_z_dot_dot(:,1,n+1) = ddx(n,3).*u1fz + ddx(n,4).*u1ez + ddx(n,5).*u2fz;

            u_y_dot(:,1,n+1) = dx(n,3).*u1fy + dx(n,4).*u1ey + dx(n,5).*u2fy;
            u_z_dot(:,1,n+1) = dx(n,3).*u1fz + dx(n,4).*u1ez + dx(n,5).*u2fz;

            u_y(:,1,n+1) = x(n,3).*u1fy + x(n,4).*u1ey + x(n,5).*u2fy;
            u_z(:,1,n+1) = x(n,3).*u1fz + x(n,4).*u1ez + x(n,5).*u2fz;

            % Transformation according to theta_p for the three modes
            theta_p_conv = deg2rad(theta_p(n+1));

            u_y_dot_dot(:,1,n+1) = u_z_dot_dot(:,1,n+1).*sin(theta_p_conv)+u_y_dot_dot(:,1,n+1).*cos(theta_p_conv);
            u_z_dot_dot(:,1,n+1) = u_z_dot_dot(:,1,n+1).*cos(theta_p_conv)-u_y_dot_dot(:,1,n+1).*sin(theta_p_conv);

            u_y_dot(:,1,n+1) = u_z_dot(:,1,n+1).*sin(theta_p_conv)+u_y_dot(:,1,n+1).*cos(theta_p_conv);
            u_z_dot(:,1,n+1) = u_z_dot(:,1,n+1).*cos(theta_p_conv)-u_y_dot(:,1,n+1).*sin(theta_p_conv);

            u_y(:,1,n+1) = u_z(:,1,n+1).*sin(theta_p_conv)+u_y(:,1,n+1).*cos(theta_p_conv);
            u_z(:,1,n+1) = u_z(:,1,n+1).*cos(theta_p_conv)-u_y(:,1,n+1).*sin(theta_p_conv);

        end

    end
    
    TI(wsp) = std(V_0_s4(5,3,:))/V_0(wsp);
    % store the converged (end) value for each V_0
    Power_converged(wsp)   = Power(end);
    omega_converged(wsp)   = omega(end);
    theta_p_converged(wsp) = theta_p(end);
    C_P_converged(wsp)     = Power_converged(wsp)./(0.5*rho*R^2*pi*V_0(wsp)^3);
    
    % store time series in difference wsp
    M_gen_Q2(:,wsp) = M_gen(1:end-1);
    Power_Q2(:,wsp) = Power(:);
    omega_Q2(:,wsp) = omega(1:N);
    theta_p_Q2(:,wsp) = theta_p(1:N);
    Thrust_Q2(:,wsp) = Thrust(:);

end

%% Rotational speed
figure(1)
plot(time,omega_Q2,'LineWidth',1)
xlabel('Time [s]','Interpreter','latex')
ylabel('Rotational speed $\omega$ [rad/s]','Interpreter','latex')
grid on;grid minor;
set(gca,'FontSize',20,'TickLabelInterpreter','latex')
set(gcf,'unit','normalized','position',[0 0 1 0.5])
saveas(gcf,'figure/omega.png')

%% Pitch angle
figure(2)
plot(time,theta_p_Q2,'LineWidth',1)
xlabel('Time [s]','Interpreter','latex')
ylabel('Pitch angle $\theta_p$ [deg]','Interpreter','latex')
grid on;grid minor;
set(gca,'FontSize',20,'TickLabelInterpreter','latex')
set(gcf,'unit','normalized','position',[0 0 1 0.5])
saveas(gcf,'figure/pitch.png')

%% Thrust
figure(3)
plot(time,Thrust_Q2,'LineWidth',1)
xlabel('Time [s]','Interpreter','latex')
ylabel('Thrust $T$ [N]','Interpreter','latex')
grid on;grid minor;
set(gca,'FontSize',20,'TickLabelInterpreter','latex')
set(gcf,'unit','normalized','position',[0 0 1 0.5])
saveas(gcf,'figure/thrust.png')

%% Blade tip deflection in flapwise and edgewise direction
figure(4)
plot(time,squeeze(u_z(end,1,1:end-1)),'LineWidth',1)
hold on
plot(time,squeeze(u_y(end,1,1:end-1)),'LineWidth',1)
hold off
ylabel('Tip Deflection [m]','Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
legend(['$u_{z,tip}$';'$u_{y,tip}$'],Interpreter='latex')
grid on; grid minor;
set(gca,'FontSize',20,'TickLabelInterpreter','latex')
set(gcf,'unit','normalized','position',[0 0 1 0.5])
saveas(gcf,'figure/tip_deflection.png')

%% PSD of Blade tip deflection in flapwise and edgewise direction

stop_lim = floor(50/delta_t);

timesim_PSD = time(stop_lim:end);
sampling_frequency = 1/delta_t;
nyquist_frequency = sampling_frequency/2;
n_parts_PSD = 2;
window = hann(round(length(timesim_PSD)/n_parts_PSD), 'periodic');
noverlap = round(length(window)/2);
nperseg = length(window);
nfft = nperseg;
% Apply Welch's method.
[sig_PSD, frequency] = pwelch(squeeze(u_z(end,1,stop_lim:end-1))-mean(squeeze(u_z(end,1,stop_lim:end-1))), window,  ...
    noverlap, nfft, sampling_frequency, 'psd');
[sig_PSD2, frequency2] = pwelch(squeeze(u_y(end,1,stop_lim:end-1))-mean(squeeze(u_y(end,1,stop_lim:end-1))), window,  ...
    noverlap, nfft, sampling_frequency, 'psd');
[sig_PSD3, frequency3] = pwelch(x_tower(stop_lim:end-1)-mean(x_tower(stop_lim:end-1)), window,  ...
    noverlap, nfft, sampling_frequency, 'psd');

figure(5)
semilogy(frequency*2*pi,sig_PSD,'LineWidth',1)
hold on
semilogy(frequency2*2*pi,sig_PSD2,'LineWidth',1)
semilogy(frequency3*2*pi,sig_PSD3,'LineWidth',1)
hold off
label = {'$\omega_{1f}$','$\omega_{1e}$','$\omega_{2f}$','$1P$','$2P$','$3P$'};
xline([omega1f omega1e omega2f omega_initial omega_initial*2 omega_initial*3] ...
    ,'--k',label,'Interpreter','latex','LineWidth',1,'LabelOrientation','horizontal','FontSize',25);
xline([1.72 7.12],'--r',{'$\omega_{tow,coup}$','$\omega_{1e,coup}$'},'Interpreter','latex','LineWidth',1,'LabelOrientation','horizontal','FontSize',25,'LabelVerticalAlignment','bottom' )
legend({'PSD($u_{z,tip}$)';'PSD($u_{y,tip}$)';'PSD($x_{tow}$)'},location="Northeast",Interpreter='latex')
xlabel('$\omega$ [rad/s]','Interpreter','latex')
xlim([0 15])
ylabel('PSD','Interpreter','latex')
grid on; grid minor;
set(gca,'FontSize',20,'TickLabelInterpreter','latex')
set(gcf,'unit','normalized','position',[0 0 1 0.5])
saveas(gcf,'figure/psd.png')
