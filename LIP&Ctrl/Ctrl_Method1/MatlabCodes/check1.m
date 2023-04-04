clear all; clc; close all;
N = 6;
is_left = false;

Lp = 0.2; 
L_min0 = -0.5;
W_min0 = 0.1; 
L_max0 = 0.5;
W_max0 = 0.4;
T_min = 0.3;
T_max = 1;
Vx = 0.8;
Vy = 0.0;

msup = 3;
mswg = 3;
mpend = 54;
mfeet = 6;
m = 60;
g = 9.8;
swingHeight = 0.1;
delta_z_vrp = 0.8;
omega = sqrt(g/delta_z_vrp);

[Tnom,Lnom,Wnom,tau_nom] = Param_fcn(L_min0, L_max0, W_min0, W_max0, T_min, T_max, omega, Vx, Vy);

% time
global t_sample
t_sample = 0.001;
t = t_sample; t0 = 0; T = Tnom; Ts = 0;

%%
if is_left
    foot_plants = [0 -(Lp/2 + Wnom) 0];
    for i=2:N+3
        foot_plants(i, :) = [(Lnom)*(i-1) (-1)^(i)*(Lp/2 + Wnom) 0];
    end
else
    foot_plants = [0 (Lp/2 + Wnom) 0];
    for i=2:N+3
        foot_plants(i, :) = [(Lnom)*(i-1) (-1)^(i-1)*(Lp/2 + Wnom) 0];
    end
end
foot_plants(end+1, :) = [1,-1, 1] .* foot_plants(end,:);

r_vrp = foot_plants;
r_vrp(:,3) =+ delta_z_vrp;

%%
% 3 Mass path generation
[r_vrp_, zmp_pend, xi_ini, xi_eos] = input3Mass(is_left, Lp, Wnom, N, Lnom, delta_z_vrp, swingHeight, Tnom, T_max, t_sample, mswg, msup, mpend, mfeet, m);
for ith = 1:N+2
    b_nom(ith,:) = (xi_eos(ith,:) - zmp_pend(ith*int32(Tnom/t_sample)+1,:));
end

% 1 Mass path generation
% [xi_ini, xi_eos] = Xi(N, r_vrp, omega, Tnom);
% for ith = 1:N+2
%     b_nom(ith,:) = (xi_eos(ith,:) - r_vrp(ith+1,:));
% end
%% initial values
u0 = [0 Lp/2]';
u0x = u0(1); u0y = u0(2);
% 3Mass IP initial pos & vel
x0_3Mass = [xi_ini(1,1) xi_ini(1,2)]';
V0_3Mass = [0 0]';
% sup leg pos array
U0_x = [];
U0_y = [];
% Swg leg final destination pos array 
UT_x = [];
UT_y = [];
SWG_traj = [];
M_FEET = [];
ZMP_FEET = [];
ZMP_PEND = [];
% CoM pos array
CoMx = [];
CoMy = [];
CoMx_3Mass = [];
CoMy_3Mass = [];
% reference DCM array
XI_ref_X = [];
XI_ref_Y = [];
XI_ref_X_ = [];
XI_ref_Y_ = [];
% measured DCM array
ZETA_mea_x = [];
ZETA_mea_y = [];
ZETA_mea_x_3Mass = [];
ZETA_mea_y_3Mass = [];
% DCM error array
ZETA_err_x = [];
ZETA_err_y = [];

PcZMP_Y = [];
PcZMP_X = [];

ini_org(1) = int32(1); fnl_org(1) = int32(Tnom/t_sample);
for i = 2:N+3
    ini_org(i) = fnl_org(i-1) + 1;
    fnl_org(i) = i*int32(Tnom/t_sample);
    if i==N+3
       fnl_org(i) =  i*int32(Tnom/t_sample)+(int32(T_max-Tnom)/t_sample);
    end
end
ini = ini_org;
fnl = fnl_org;
s = 1; inc(s) = 0;

Step = 1; i = 1; q = 1; n = 0; 
qpresult = [0;0;0;0;0];
F = 0; 
%% control loop
while Step(i) == 1
    s = s + 1;

    % Disturbance insertation
    if n+1 == 3 && t <= 0.1
        F = 120; %Max 290 for 1 Mass path generaton
    else
        F = 0;
    end

    %% regenerate DCM pattern 
    % 3Mass path generation
    xi_X = zmp_pend(ini(n+1)-1+q,1) + exp(omega*(t-T))*(zmp_pend(ini(n+2),1) + b_nom(n+1,1) - zmp_pend(ini(n+1)-1+q,1));
    xi_Y = zmp_pend(ini(n+1)-1+q,2) + exp(omega*(t-T))*(zmp_pend(ini(n+2),2) + b_nom(n+1,2) - zmp_pend(ini(n+1)-1+q,2));
    % 1Mass path generation
%     xi_X = r_vrp(n+1,1) + exp(omega*(t-T))*(r_vrp(n+2,1) + b_nom(n+1,1) - r_vrp(n+1,1));
%     xi_Y = r_vrp(n+1,2) + exp(omega*(t-T))*(r_vrp(n+2,2) + b_nom(n+1,2) - r_vrp(n+1,2));

    XI_ref_X = horzcat(XI_ref_X,[t + sum(Ts) xi_X]');
    XI_ref_Y = horzcat(XI_ref_Y,[t + sum(Ts) xi_Y]');
    
    %% measured com and dcm of 3Mass IP

    k1x = t_sample*f1(t0,x0_3Mass(1),V0_3Mass(1));
    l1x = t_sample*f2(t0,x0_3Mass(1),V0_3Mass(1), zmp_pend(ini(n+1)-1+q,1),0);

    k2x = t_sample*f1(t0+t_sample/2, x0_3Mass(1)+k1x/2, V0_3Mass(1)+l1x/2);
    l2x = t_sample*f2(t0+t_sample/2, x0_3Mass(1)+k1x/2, V0_3Mass(1)+l1x/2, zmp_pend(ini(n+1)-1+q,1),0);

    k3x = t_sample*f1(t0+t_sample/2, x0_3Mass(1)+k2x/2, V0_3Mass(1)+l2x/2);
    l3x = t_sample*f2(t0+t_sample/2, x0_3Mass(1)+k2x/2, V0_3Mass(1)+l2x/2, zmp_pend(ini(n+1)-1+q,1),0);

    k4x = t_sample*f1(t0+t_sample, x0_3Mass(1)+k3x, V0_3Mass(1)+l3x);
    l4x = t_sample*f2(t0+t_sample, x0_3Mass(1)+k3x, V0_3Mass(1)+l3x, zmp_pend(ini(n+1)-1+q,1),0);
    
    k1y = t_sample*f1(t0,x0_3Mass(2),V0_3Mass(2));
    l1y = t_sample*f2(t0,x0_3Mass(2),V0_3Mass(2), zmp_pend(ini(n+1)-1+q,2),F);

    k2y = t_sample*f1(t0+t_sample/2, x0_3Mass(2)+k1y/2, V0_3Mass(2)+l1y/2);
    l2y = t_sample*f2(t0+t_sample/2, x0_3Mass(2)+k1y/2, V0_3Mass(2)+l1y/2, zmp_pend(ini(n+1)-1+q,2),F);

    k3y = t_sample*f1(t0+t_sample/2, x0_3Mass(2)+k2y/2, V0_3Mass(2)+l2y/2);
    l3y = t_sample*f2(t0+t_sample/2, x0_3Mass(2)+k2y/2, V0_3Mass(2)+l2y/2, zmp_pend(ini(n+1)-1+q,2),F);

    k4y = t_sample*f1(t0+t_sample, x0_3Mass(2)+k3y, V0_3Mass(2)+l3y);
    l4y = t_sample*f2(t0+t_sample, x0_3Mass(2)+k3y, V0_3Mass(2)+l3y, zmp_pend(ini(n+1)-1+q,2),F);
    
    % update values for t0
    V0_3Mass = V0_3Mass + t_sample*(omega^2*(x0_3Mass-zmp_pend(ini(n+1)-1+q,1:2)')+[0;-F/m]);
    x0_3Mass(1) = x0_3Mass(1) + (k1x + 2*k2x + 2*k3x + k4x)/6;
    x0_3Mass(2) = x0_3Mass(2) + (k1y + 2*k2y + 2*k3y + k4y)/6;
    t0 = t0 + t_sample;

    % calculate com and dcm for t0
    CoM_x_3Mass = [t0 + sum(Ts) x0_3Mass(1)]';
    CoM_y_3Mass = [t0 + sum(Ts) x0_3Mass(2)]';
    CoMx_3Mass = horzcat(CoMx_3Mass,CoM_x_3Mass);
    CoMy_3Mass = horzcat(CoMy_3Mass,CoM_y_3Mass);
    
    xi_meas_3Mass = x0_3Mass + V0_3Mass/omega;
    zeta_mea_x = [t + sum(Ts) xi_meas_3Mass(1)]'; %simoutx(q,2) xi_meas_3Mass(1)
    zeta_mea_y = [t + sum(Ts) xi_meas_3Mass(2)]'; %simouty(q,2) xi_meas_3Mass(2)
    ZETA_mea_x_3Mass = horzcat(ZETA_mea_x_3Mass, zeta_mea_x);
    ZETA_mea_y_3Mass = horzcat(ZETA_mea_y_3Mass, zeta_mea_y);
    
    %% dcm error 
    zeta_err_x = [t + sum(Ts) xi_meas_3Mass(1)-xi_X]'; %simoutx(q,2) xi_meas_3Mass(1)
    zeta_err_y = [t + sum(Ts) xi_meas_3Mass(2)-xi_Y]'; %simouty(q,2) xi_meas_3Mass(2)
    ZETA_err_x = horzcat(ZETA_err_x,zeta_err_x);
    ZETA_err_y = horzcat(ZETA_err_y,zeta_err_y);
    
    %% QP 
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
    PcZMP_Y = horzcat(PcZMP_Y, PcZMP_y(q,n+1)+u0y);
    PcZMP_X = horzcat(PcZMP_X, PcZMP_x(q,n+1)+u0x);
    
    %% update QP constraint parameters
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
          zeta_err_x, zeta_err_y, PcZMP_y(q,n+1), PcZMP_x(q,n+1)); % zeta_err_x, zeta_err_y, PcZMP_y(q,n+1), PcZMP_x(q,n+1) [0 0], [0 0], 0, 0

    T = (1/omega)*log(Opt_Vector(3));
    
    % Swg leg new destination
    uT_x = [t + sum(Ts) Opt_Vector(1)]';
    uT_y = [t + sum(Ts) Opt_Vector(2)]';
    UT_x = horzcat(UT_x, uT_x);
    UT_y = horzcat(UT_y, uT_y);
    
    % Sup leg pos
    u0_x = [t + sum(Ts) u0x]';
    u0_y = [t + sum(Ts) u0y]';
    U0_x = horzcat(U0_x, u0_x);
    U0_y = horzcat(U0_y, u0_y);
    
    %% update pattern parameter
    for j = n:N+2
       r_vrp(j+2,1) = r_vrp(j+2,1) + qpresult(1);
       r_vrp(j+2,2) = r_vrp(j+2,2) + qpresult(3);
    end
    
    if is_left
        r_f_l = [1 1 0].*r_vrp;
        r_f_r = [1 1 0].*r_vrp;
        r_f_l(1,:) = [0 (Lp/2 + Wnom) 0];
        r_f_l(3:2:N+3,:)=r_f_l(2:2:N+2,:);
        r_f_r(2:2:N+3,:)=r_f_r(1:2:N+2,:);
    else
        r_f_r = [1 1 0].*r_vrp;
        r_f_l = [1 1 0].*r_vrp;
        r_f_r(1,:) = [0 -(Lp/2 + Wnom) 0];
        r_f_r(3:2:N+3,:)=r_f_r(2:2:N+2,:);
        r_f_l(2:2:N+3,:)=r_f_l(1:2:N+2,:);
    end
    inc(s) = floor(T/t_sample) - int32(Tnom/t_sample);
    for ith = n+1:N+3
        ini(ith+1:end) = ini_org(ith+1:end) + int32(inc(s)*ones(1,N+3-ith));
        fnl(ith:end) = fnl_org(ith:end) + int32(inc(s)*ones(1,N+4-ith));
    end  
    %% going next step
    t = t + t_sample;
    
    [Opt_Vector(1); Opt_Vector(2); Opt_Vector(3); Opt_Vector(4); Opt_Vector(5)]
    [n T t q]
    
    if t>T
        t = t_sample;
        t0 = 0;
        Ts(i) = T;
        i = i+1;
        Step(i) = 1;
        n = length(Step)-1;
        
        % new Sup leg pos = last Swg leg destination pos
        u0x = Opt_Vector(1);
        u0y = Opt_Vector(2);
        u0 = [u0x u0y]';
        q = 0;
        
        ini_org = ini;
        fnl_org = fnl;
    end
    
    %% check ending condition
    if n == N+2 % N
       break
    end
    
    q = q+1;
    
    
    %% swg leg traj generation
    for ith = n+1:n+2 % calculate zmp_pend for current step and next step
        init_time = double(ini(ith))*t_sample;
        final_time = double(fnl(ith))*t_sample;
        if is_left
            if mod(ith,2) ~= 0
                swingfootpos0 = r_f_l(ith, :);
                swingfootpos1 = r_f_l(ith+1, :);
            else
                swingfootpos0 = r_f_r(ith, :);
                swingfootpos1 = r_f_r(ith+1, :);
            end
        else
            if mod(ith,2) ~= 0
                swingfootpos0 = r_f_r(ith, :);
                swingfootpos1 = r_f_r(ith+1, :);
            else
                swingfootpos0 = r_f_l(ith, :);
                swingfootpos1 = r_f_l(ith+1, :);
            end
        end

        [qswing, dqswing, ddqswing] = getSwingFootTraj(swingfootpos0', swingfootpos1', swingHeight, ...
                            init_time, final_time,t_sample);   

        swg_traj(ini(ith):fnl(ith),:) = [qswing' dqswing' ddqswing'];
%         SWG_traj = horzcat(SWG_traj,swg_traj);

        %% three mass parameters evaluation 
        m_feet(ini(ith):fnl(ith),:) = mswg*(qswing(1:2,:)' -  ones(length(qswing),1)*[r_vrp(ith,1) r_vrp(ith,2)]).*(g*ones(length(qswing),1) + ddqswing(3,:)')...
            - mswg*(ddqswing(1:2,:)').*(qswing(3,:)' - delta_z_vrp*ones(length(qswing),1));
%         M_FEET = horzcat(M_FEET, m_feet);

        zmp_feet(ini(ith):fnl(ith),:) = m_feet(ini(ith):fnl(ith),:)/(mfeet*g) + ones(length(qswing),1)*[r_vrp(ith,1) r_vrp(ith,2)];
%         ZMP_FEET = horzcat(ZMP_FEET, zmp_feet);

        zmp_pend(ini(ith):fnl(ith),1:2) = (m/mpend)*ones(length(qswing),1)*[r_vrp(ith,1) r_vrp(ith,2)] - (mfeet/mpend)*zmp_feet(ini(ith):fnl(ith),:);

    end
    ZMP_PEND = horzcat(ZMP_PEND, [t + sum(Ts); zmp_pend(ini(n+1)-1+q,1:2)']);
end
%% plot result
figure(1)
plot(XI_ref_X(1,:),XI_ref_X(2,:),'color','k','LineStyle','--','linewidth',1.5);hold on;
plot(ZETA_mea_x_3Mass(1,:),ZETA_mea_x_3Mass(2,:),'color','k','linewidth',1.5);hold on;
plot(CoMx_3Mass(1,:),CoMx_3Mass(2,:),'color','m','linewidth',1.5);hold on;
plot(UT_x(1,:),UT_x(2,:),'color','b','linewidth',1.5);hold on;
plot(U0_x(1,:),U0_x(2,:),'color','c','linewidth',1.5);hold on;
plot(ZMP_PEND(1,:), ZMP_PEND(2,:),'color','g','LineStyle','-','linewidth',1.5)
legend('\xi_{ref,x}','\xi_{meas,x}','x_{com,meas}','u_{T,x}','u_{0,x}','zmp_{pend}') %'P_{cZMP,x}','\xi_{meas,x}','x_{com,meas}','u_{T,x}','swg_{x}','swg_{z}','zmp_{feet}',
%%
figure(2)
plot(XI_ref_Y(1,:),XI_ref_Y(2,:),'color','r','linewidth',1.5);hold on;
plot(ZETA_mea_y_3Mass(1,:),ZETA_mea_y_3Mass(2,:),'color','k','linewidth',1.5);hold on;
plot(CoMy_3Mass(1,:),CoMy_3Mass(2,:),'color','m','linewidth',1.5);hold on;
plot(UT_y(1,:),UT_y(2,:),'color','b','linewidth',1.5);hold on;
plot(U0_y(1,:),U0_y(2,:),'color','c','linewidth',1.5);
plot(ZMP_PEND(1,:), ZMP_PEND(3,:),'color','g','LineStyle','-','linewidth',1.5)
legend('\xi_{ref,y}','\xi_{meas,y}','y_{com,meas}','u_{T,y}','u_{0,y}','zmp_{pend}') %,'\xi_{meas,y}','y_{com,meas}','u_{T,y}','swg_{y}','zmp_{feet}','P_{cZMP,y}'

%functions definition
function [xi_ini, xi_eos] = Xi(N, r_vrp, omega, Tnom)
xi_eos = zeros(N+2,3);
xi_eos(N+2,:) = r_vrp(end-1,:);
for i = N+1:-1:1
        xi_eos(i,:) = r_vrp(i+1,:) + (exp(-omega*Tnom))*(xi_eos(i+1,:)-r_vrp(i+1,:));
        xi_ini(i+1,:) = xi_eos(i,:);
end
xi_ini(1,:) = r_vrp(i,:) + (exp(-omega*Tnom))*(xi_eos(i,:)-r_vrp(i,:));
end

function dxdt = f1(t,x,v)
dxdt = v;
end

function dvdt = f2(t,x,v,u,F)
omega = 3.5;
m = 60;
d = -F/m;
dvdt = omega^2*(x-u) + d;
end