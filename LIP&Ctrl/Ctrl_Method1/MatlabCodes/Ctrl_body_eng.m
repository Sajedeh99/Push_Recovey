clear all; clc; close all;
N = 2;
is_left = false;

Lp = 0.2;
L_min0 = -0.5;
W_min0 = 0.1; 
L_max0 = 0.5;
W_max0 = 0.4;
T_min = 0.3;
T_max = 1;
Vx = 0.5;
Vy = 0.0;

m = 60;
g = 9.8; 
delta_z_vrp = 0.8;
omega = sqrt(g/delta_z_vrp);

[Tnom,Lnom,Wnom,tau_nom] = Param_fcn(L_min0, L_max0, W_min0, W_max0, T_min, T_max, omega, Vx, Vy);
global t_sample
t_sample = 0.001;

traj_ds = struct('index', {}, 'xi_ds', {}, 'state', {});
traj_ss = struct('index', {}, 'xi_ss', {}, 'state', {});

% for i=1:N+3
%     if i == 1
%         r_f_r(i,1) = 0;
%         r_f_l(i,1) = 0;
%     elseif i == N+3
%         r_f_l(i, 1) = (Lnom)*(i-2);
%         r_f_r(i,1) = (Lnom)*(i-2);
%     elseif mod(i,2) == 0
%         r_f_l(i,1) = (Lnom)*(i-1);
%         r_f_r(i,1) = r_f_r(i-1,1);  
%     else
%         r_f_r(i,1) = (Lnom)*(i-1);
%         r_f_l(i,1) = r_f_l(i-1,1);
%     end
%     r_f_l(i,2:3) = [(Lp/2 + Wnom) 0];
%     r_f_r(i,2:3) = [-(Lp/2 + Wnom) 0];
% end
% r_vrp = r_f_r;
% r_vrp(2:2:N+3, 1:2) = r_f_l(2:2:N+3,1:2);
% if is_left == false
%     temp = r_f_l;
%     r_f_l = [1 -1 1].*r_f_r;
%     r_f_r = [1 -1 1].*temp;
%     r_vrp = r_f_l;
%     r_vrp(2:2:N+3, 1:2) = r_f_r(2:2:N+3,1:2);
% end
% r_vrp(:,3) =+ delta_z_vrp;
% r_vrp(end,:) = [1 1 1].*r_vrp(end,:);

%%
if is_left
    foot_plants = [0 -(Lp/2 + Wnom) 0];
    for i=2:N+2
        foot_plants(i, :) = [(Lnom)*(i-1) (-1)^(i)*(Lp/2 + Wnom) 0];
    end
else
    foot_plants = [0 (Lp/2 + Wnom) 0];
    for i=2:N+2
        foot_plants(i, :) = [(Lnom)*(i-1) (-1)^(i-1)*(Lp/2 + Wnom) 0];
    end
end
foot_plants(end+1, :) = [1,-1, 1] .* foot_plants(end,:);

r_vrp = foot_plants;
r_vrp(:,3) =+ delta_z_vrp;

%%

[xi_ini, xi_eos] = Xi(N, r_vrp, omega, Tnom);
for ith = 1:N+2
    b_nom(ith,:) = (xi_eos(ith,:) - r_vrp(ith+1,:));
end
%% initial values
% Sup leg initial pos
u0 = [0 Lp/2]';
u0x = u0(1);
u0y = u0(2);
% IP initial pos & vel
x0 = [xi_ini(1,1) xi_ini(1,2)]';
V0 = [0 0]';
% sup leg pos array
U0_x = [];
U0_y = [];
% Swg leg final destination pos array 
UT_x = [];
UT_y = [];
% CoM pos array
CoMx = [];
CoMy = [];
PcZMP_Y = [];
PcZMP_X = [];
% reference DCM array
XI_ref_X = [];
XI_ref_Y = [];
% measured DCM array
ZETA_mea_x = [];
ZETA_mea_y = [];
% DCM error array
ZETA_err_x = [];
ZETA_err_y = [];
% time 
t = t_sample; T = Tnom; Ts = 0; Tsim = [t_sample:t_sample:T_max];
Step = 1; i = 1; q = 0; n = 0; s = 0;
qpresult = [0;0;0;0;0];
init_time = t
final_time = T
%% control loop
while Step(i) == 1
    q = q+1;
    s = s+1;
    % update pattern parameter
    for j = n:N+1
       r_vrp(j+2,1) = r_vrp(j+2,1) + qpresult(1);
       r_vrp(j+2,2) = r_vrp(j+2,2) + qpresult(3);
    end
    
    if is_left
        r_f_l = [1 1 0].*r_vrp;
        r_f_r = [1 1 0].*r_vrp;
        r_f_l(1,:) = [0 (Lp/2 + Wnom) 0];
        r_f_l(3:2:N+3,:)=r_f_l(2:2:N+floor((N+2)/2),:);
        r_f_r(2:2:N+3,:)=r_f_r(1:2:N+floor((N+2)/2),:);
    else
        r_f_r = [1 1 0].*r_vrp;
        r_f_l = [1 1 0].*r_vrp;
        r_f_r(1,:) = [0 -(Lp/2 + Wnom) 0];
        r_f_r(3:2:N+3,:)=r_f_r(2:2:N+floor((N+2)/2),:);
        r_f_l(2:2:N+3,:)=r_f_l(1:2:N+floor((N+2)/2),:);
    end
    
    % regenerate DCM pattern 
    xi_X = r_vrp(n+1,1) + exp(omega*(t-T))*(r_vrp(n+2,1) + b_nom(n+1,1) - r_vrp(n+1,1));
    xi_Y = r_vrp(n+1,2) + exp(omega*(t-T))*(r_vrp(n+2,2) + b_nom(n+1,2) - r_vrp(n+1,2));
    
    % simulate IP with initial value (u0, x0, v0) of ith Step
    if q == 1
        sim('LIPM_Dynamicsx',[t_sample T_max]);
        sim('LIPM_Dynamicsy',[t_sample T_max]);
    end
    
    % time variable is the local time
    time = Tsim(q) + sum(Ts);
    
    % measured com and dcm of IP
    CoM_x = [time simoutx(q,1)]';
    CoM_y = [time simouty(q,1)]';
    CoMx = horzcat(CoMx,CoM_x);
    CoMy = horzcat(CoMy,CoM_y);
    zeta_mea_x = [time simoutx(q,2)]';
    zeta_mea_y = [time simouty(q,2)]';
    ZETA_mea_x = horzcat(ZETA_mea_x,zeta_mea_x);
    ZETA_mea_y = horzcat(ZETA_mea_y,zeta_mea_y);
    
    % reference dcm
    xi_ref_X = [time xi_X]';
    xi_ref_Y = [time xi_Y]';
    XI_ref_X = horzcat(XI_ref_X,xi_ref_X);
    XI_ref_Y = horzcat(XI_ref_Y,xi_ref_Y);
    
    % dcm error 
    zeta_err_x = [time simoutx(q,2)-xi_X]';
    zeta_err_y = [time simouty(q,2)-xi_Y]';
    ZETA_err_x = horzcat(ZETA_err_x,zeta_err_x);
    ZETA_err_y = horzcat(ZETA_err_y,zeta_err_y);
    
    % update QP constraint parameters 
    PcZMP_y(q,n+1) = -(exp(omega*(T-t+0.05)))*zeta_err_y(2)/(1-exp(omega*(T-t+0.05)));
    PcZMP_x(q,n+1) = -(exp(omega*(T-t+0.02)))*zeta_err_x(2)/(1-exp(omega*(T-t+0.02)));
    
    if abs(PcZMP_y(q,n+1)) >= 0.04
       if  PcZMP_y(q,n+1) > 0
           PcZMP_y(q,n+1) = 0.04;
       else
           PcZMP_y(q,n+1) = -0.04;
       end
    end
    if abs(PcZMP_x(q,n+1)) >= 0.08
       if  PcZMP_x(q,n+1) > 0
           PcZMP_x(q,n+1) = 0.08;
       else
           PcZMP_x(q,n+1) = -0.08;
       end
    end 
    L_min = u0(1) - .5;
    L_max = u0(1) + .5;
    if mod(n,2) == 0
        W_min = u0(2) - W_max0;
        W_max = u0(2) - W_min0;
    else
        W_min = u0(2) + W_min0;
        W_max = u0(2) + W_max0;
    end
    
    % QP
    [qpresult, Opt_Vector] = controller_eng(t, T, Lnom, Wnom, L_min, L_max, W_min, W_max, T_min, T_max,...
          b_nom(n+1,1), b_nom(n+1,2), omega, zeta_mea_x, zeta_mea_y, r_vrp(n+2,1), r_vrp(n+2,2),...
          zeta_err_x, zeta_err_y, PcZMP_y(q,n+1), PcZMP_x(q,n+1)); %PcZMP_y(q,n+1), PcZMP_x(q,n+1)

    T = (1/omega)*log(Opt_Vector(3));
 
    % swg foot traj generation
    final_time = T + sum(Ts)
    
    % Sup leg pos
    u0_x = [t + sum(Ts) u0x]';
    u0_y = [t + sum(Ts) u0y]';
    U0_x = horzcat(U0_x, u0_x);
    U0_y = horzcat(U0_y, u0_y);
    
    % Swg leg new destination
    uT_x = [t + sum(Ts) Opt_Vector(1)]';
    uT_y = [t + sum(Ts) Opt_Vector(2)]';
    UT_x = horzcat(UT_x, uT_x);
    UT_y = horzcat(UT_y, uT_y);

    PcZMP_Y = horzcat(PcZMP_Y, PcZMP_y(q,n+1)+u0y);
    PcZMP_X = horzcat(PcZMP_X, PcZMP_x(q,n+1)+u0x);
    
    t = t + t_sample;
    
%     [Opt_Vector(1); Opt_Vector(2); Opt_Vector(3); Opt_Vector(4); Opt_Vector(5)]
    [n T t]
    
    % going next step
    if t>T
        t = 0;
        Ts(i) = T;
        init_time = sum(Ts)
        i = i+1;
        Step(i) = 1;
        n = length(Step)-1;
        % initial x0 & v0 for IP in next step
        x0 = [CoM_x(2) CoM_y(2)]';
        term1 = [simoutx(q,2), simouty(q,2)];
        term2 = [simoutx(q,1), simouty(q,1)];
        V0 = omega*(term1-term2);
        
        % new Sup leg pos = last Swg leg destination pos
        u0x = Opt_Vector(1);
        u0y = Opt_Vector(2);
        u0 = [u0x u0y]';
        q = 0;
    end

    if n == N+2 % N
        break
    end
end
%% plot result
figure(1)
plot(XI_ref_X(1,:),XI_ref_X(2,:),'color','k','LineStyle','-','linewidth',2);hold on;
plot(ZETA_mea_x(1,:),ZETA_mea_x(2,:),'color','g');hold on;
plot(CoMx(1,:),CoMx(2,:),'color','m');hold on;
plot(UT_x(1,:),UT_x(2,:),'color','b');hold on;
plot(U0_x(1,:),U0_x(2,:),'color','c','linewidth',2);
% plot(ZETA_mea_x(1,:),PcZMP_X,'color','r','linewidth',2);
legend('\xi_{ref,x}','\xi_{meas,x}','x_{com,meas}','u_{T,x}','u_{0,x}') %,'P_{cZMP,x}'
figure(2)
plot(XI_ref_Y(1,:),XI_ref_Y(2,:),'color','k','LineStyle','-','linewidth',2);hold on;
plot(ZETA_mea_y(1,:),ZETA_mea_y(2,:),'color','g','linewidth',2);hold on;
plot(CoMy(1,:),CoMy(2,:),'color','m','linewidth',2);hold on;
plot(UT_y(1,:),UT_y(2,:),'color','b','linewidth',2);hold on;
plot(U0_y(1,:),U0_y(2,:),'color','c','linewidth',2);
% plot(ZETA_mea_y(1,:),PcZMP_Y,'color','r','linewidth',2);
legend('\xi_{ref,y}','\xi_{meas,y}','y_{com,meas}','u_{T,y}','u_{0,y}') %,'P_{cZMP,y}'

%functions definition
function [xi_ini, xi_eos] = Xi(N, r_vrp, omega, Tnom)
xi_eos = zeros(N+2,3);
xi_eos(N+2,:) = r_vrp(end,:);
for i = N+1:-1:1
        xi_eos(i,:) = r_vrp(i+1,:) + (exp(-omega*Tnom))*(xi_eos(i+1,:)-r_vrp(i+1,:));
        xi_ini(i+1,:) = xi_eos(i,:);
end
xi_ini(1,:) = r_vrp(i,:) + (exp(-omega*Tnom))*(xi_eos(i,:)-r_vrp(i,:));
end
