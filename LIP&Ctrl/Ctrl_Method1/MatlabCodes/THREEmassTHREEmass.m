clear all; clc; close all;
N = 9;
is_left = false;

animateOn = true; 
speedupfactor = 15; % animation speed up 

Lp = 0.2; 
L_min0 = -0.5;
W_min0 = 0.1; 
L_max0 = 0.5;
W_max0 = 0.4;
T_min = 0.3;
T_max = 1;
Vx = 1.5;
Vy = 0.0;

msup = 3;
mswg = 3;
mpend = 54;
mfeet = 6;
m = 60;
g = 9.8;
swingHeight = 0.1;
z_robot = 0.68;
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
foot_plants(end+1, :) = [1,-0, 1] .* foot_plants(end,:);

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
xi_meas_3Mass = x0_3Mass + V0_3Mass/omega;
% sup leg pos array
U0_x = [];
U0_y = [];
% Swg leg final destination pos array 
UT_x = []; UT_xEr = [];
UT_y = []; UT_yEr = [];
bT_xEr = [];
bT_yEr = [];
T_step = []; t_var = [];
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

PcZMP_Y = []; PcZMP_XX = [];
PcZMP_X = []; PcZMP_YY = [];

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
xi_X = r_vrp(n+1,1) + exp(omega*(t-T))*(r_vrp(n+2,1) + b_nom(n+1,1) - r_vrp(n+1,1));
xi_Y = r_vrp(n+1,2) + exp(omega*(t-T))*(r_vrp(n+2,2) + b_nom(n+1,2) - r_vrp(n+1,2));

x_com(1,:) = [xi_ini(1,1) xi_ini(1,2)];
com_dot = [0, 0];
stateR(:,s) = [-u0y; u0x; -0];
stateL(:,s) = [ u0y; u0x; -0];
stateCoM(:,s)=[x0_3Mass(2); x0_3Mass(1); z_robot];
% [robot, hLeftRel, hRightRel, hCoMRel] = createRobot(x0_3Mass, z_robot, stateL, stateR, stateCoM);
% animate(stateR(:,s), stateL(:,s), stateCoM(:,s), animateOn, robot, hLeftRel, hRightRel, hCoMRel, speedupfactor, s, z_robot);

%% control loop
while Step(i) == 1
    s = s + 1;
% [xi_meas_3Mass(2); xi_Y]
    % Disturbance insertation
    if n+1 == 3 && t <= 0.1
%         if t<=0.1 & t>= 0.09
%             plot3([x0_3Mass(2) x0_3Mass(2)+0.25],[x0_3Mass(1) x0_3Mass(1)],[z_robot+0.2 z_robot+0.2],'LineWidth',2.5,'Color','m');
%             plot3([x0_3Mass(2) x0_3Mass(2)+0.05],[x0_3Mass(1) x0_3Mass(1)+0.05],[z_robot+0.2 z_robot+0.2],'LineWidth',2.5,'Color','m');
%             plot3([x0_3Mass(2) x0_3Mass(2)+0.05],[x0_3Mass(1) x0_3Mass(1)-0.05],[z_robot+0.2 z_robot+0.2],'LineWidth',2.5,'Color','m');
%         end
        Fy = 360; % max 380 @ 0.1s
        Fx = 0;
    else
        Fy = 0;
        Fx = 0;
    end
    
    %% reference dcm and com
    XI_ref_X = horzcat(XI_ref_X,[t + sum(Ts) xi_X]');
    XI_ref_Y = horzcat(XI_ref_Y,[t + sum(Ts) xi_Y]');

    x_com(s,:) = x_com(s-1,:) - t_sample*omega*(x_com(s-1,:) - [xi_X xi_Y]);

    %% regenerate DCM pattern 
    % 3Mass path generation
    xi_X = zmp_pend(ini(n+1)-1+q,1) + exp(omega*(t-T))*(zmp_pend(ini(n+2),1) + b_nom(n+1,1) - zmp_pend(ini(n+1)-1+q,1));
    xi_Y = zmp_pend(ini(n+1)-1+q,2) + exp(omega*(t-T))*(zmp_pend(ini(n+2),2) + b_nom(n+1,2) - zmp_pend(ini(n+1)-1+q,2));
    % 1Mass path generation
%     xi_X = r_vrp(n+1,1) + exp(omega*(t-T))*(r_vrp(n+2,1) + b_nom(n+1,1) - r_vrp(n+1,1));
%     xi_Y = r_vrp(n+1,2) + exp(omega*(t-T))*(r_vrp(n+2,2) + b_nom(n+1,2) - r_vrp(n+1,2));

    com_dot(s, :) = (omega) * ([xi_X xi_Y] - x_com(s,:)); %  (omega^2) was wrong
    
    %% dcm error 
    zeta_err_x = [t + sum(Ts) xi_meas_3Mass(1)-xi_X]'; %simoutx(q,2) xi_meas_3Mass(1)
    zeta_err_y = [t + sum(Ts) xi_meas_3Mass(2)-xi_Y]'; %simouty(q,2) xi_meas_3Mass(2)
    ZETA_err_x = horzcat(ZETA_err_x,zeta_err_x);
    ZETA_err_y = horzcat(ZETA_err_y,zeta_err_y);
    
    % DCM end-of-step controller
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
    PcZMP_YY = horzcat(PcZMP_YY, PcZMP_y(q,n+1));
    PcZMP_XX = horzcat(PcZMP_XX, PcZMP_x(q,n+1));
     
    %% measured com and dcm of 3Mass IP
    k1x = t_sample*f1(t0,x0_3Mass(1),V0_3Mass(1));
    l1x = t_sample*f2(t0,x0_3Mass(1),V0_3Mass(1), zmp_pend(ini(n+1)-1+q,1)+PcZMP_x(q,n+1),Fx);

    k2x = t_sample*f1(t0+t_sample/2, x0_3Mass(1)+k1x/2, V0_3Mass(1)+l1x/2);
    l2x = t_sample*f2(t0+t_sample/2, x0_3Mass(1)+k1x/2, V0_3Mass(1)+l1x/2, zmp_pend(ini(n+1)-1+q,1)+PcZMP_x(q,n+1),Fx);

    k3x = t_sample*f1(t0+t_sample/2, x0_3Mass(1)+k2x/2, V0_3Mass(1)+l2x/2);
    l3x = t_sample*f2(t0+t_sample/2, x0_3Mass(1)+k2x/2, V0_3Mass(1)+l2x/2, zmp_pend(ini(n+1)-1+q,1)+PcZMP_x(q,n+1),Fx);

    k4x = t_sample*f1(t0+t_sample, x0_3Mass(1)+k3x, V0_3Mass(1)+l3x);
    l4x = t_sample*f2(t0+t_sample, x0_3Mass(1)+k3x, V0_3Mass(1)+l3x, zmp_pend(ini(n+1)-1+q,1)+PcZMP_x(q,n+1),Fx);
    
    k1y = t_sample*f1(t0,x0_3Mass(2),V0_3Mass(2));
    l1y = t_sample*f2(t0,x0_3Mass(2),V0_3Mass(2), zmp_pend(ini(n+1)-1+q,2)+PcZMP_y(q,n+1),Fy);

    k2y = t_sample*f1(t0+t_sample/2, x0_3Mass(2)+k1y/2, V0_3Mass(2)+l1y/2);
    l2y = t_sample*f2(t0+t_sample/2, x0_3Mass(2)+k1y/2, V0_3Mass(2)+l1y/2, zmp_pend(ini(n+1)-1+q,2)+PcZMP_y(q,n+1),Fy);

    k3y = t_sample*f1(t0+t_sample/2, x0_3Mass(2)+k2y/2, V0_3Mass(2)+l2y/2);
    l3y = t_sample*f2(t0+t_sample/2, x0_3Mass(2)+k2y/2, V0_3Mass(2)+l2y/2, zmp_pend(ini(n+1)-1+q,2)+PcZMP_y(q,n+1),Fy);

    k4y = t_sample*f1(t0+t_sample, x0_3Mass(2)+k3y, V0_3Mass(2)+l3y);
    l4y = t_sample*f2(t0+t_sample, x0_3Mass(2)+k3y, V0_3Mass(2)+l3y, zmp_pend(ini(n+1)-1+q,2)+PcZMP_y(q,n+1),Fy);
    
    % update values for t0
    V0_3Mass = V0_3Mass + t_sample*(omega^2*(x0_3Mass-(zmp_pend(ini(n+1)-1+q,1:2)'+[PcZMP_x(q,n+1);PcZMP_y(q,n+1)]))+[-Fx/m;-Fy/m]);
    x0_3Mass(1) = x0_3Mass(1) + (k1x + 2*k2x + 2*k3x + k4x)/6;
    x0_3Mass(2) = x0_3Mass(2) + (k1y + 2*k2y + 2*k3y + k4y)/6;
    t0 = t0 + t_sample;

    % calculate com and dcm for t0
    CoM_x_3Mass = [t0 + sum(Ts) x0_3Mass(1)]';
    CoM_y_3Mass = [t0 + sum(Ts) x0_3Mass(2)]';
    CoMx_3Mass = horzcat(CoMx_3Mass,CoM_x_3Mass);
    CoMy_3Mass = horzcat(CoMy_3Mass,CoM_y_3Mass);
    
    stateCoM(:,s) = [x0_3Mass(2);x0_3Mass(1);z_robot];

    % measured DCM pattern
    xi_meas_3Mass = x0_3Mass + V0_3Mass/omega;
    zeta_mea_x = [t + sum(Ts) xi_meas_3Mass(1)]'; %simoutx(q,2) xi_meas_3Mass(1)
    zeta_mea_y = [t + sum(Ts) xi_meas_3Mass(2)]'; %simouty(q,2) xi_meas_3Mass(2)
    ZETA_mea_x_3Mass = horzcat(ZETA_mea_x_3Mass, zeta_mea_x);
    ZETA_mea_y_3Mass = horzcat(ZETA_mea_y_3Mass, zeta_mea_y);

    
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
% [xi_meas_3Mass(2); xi_Y]
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
    
    UT_xEr = horzcat(UT_xEr, [t + sum(Ts) Opt_Vector(1)-foot_plants(n+2,1)]');
    UT_yEr = horzcat(UT_yEr, [t + sum(Ts) Opt_Vector(2)-foot_plants(n+2,2)]');
    bT_xEr = horzcat(bT_xEr, [t + sum(Ts) qpresult(2)]');
    bT_yEr = horzcat(bT_yEr, [t + sum(Ts) qpresult(4)]');
    T_step = horzcat(T_step, [t + sum(Ts) T]');
    t_var = horzcat(t_var, [t + sum(Ts) t]');
    
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
    [n T t]
    
    if t>T
        t = t_sample;
        t0 = 0;
        Ts(i) = T;
        i = i+1;
        Step(i) = 1;
        n = length(Step)-1;
        
        % new Sup leg pos = last Swg leg destination pos
        foot_plants = r_vrp;
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
    if (is_left==false)
        if mod(n,2)==0
            stateR(:,s) = [swg_traj(ini(n+1)-1+q,2); swg_traj(ini(n+1)-1+q,1); swg_traj(ini(n+1)-1+q,3)];
            stateL(:,s) = [u0y; u0x; -0];
        else
            stateR(:,s) = [u0y; u0x; -0];
            stateL(:,s) = [swg_traj(ini(n+1)-1+q,2); swg_traj(ini(n+1)-1+q,1); swg_traj(ini(n+1)-1+q,3)];
        end        
    else
        if mod(n,2)==0
            stateR(:,s) = [u0y; u0x; -0];
            stateL(:,s) = [swg_traj(ini(n+1)-1+q,2); swg_traj(ini(n+1)-1+q,1); swg_traj(ini(n+1)-1+q,3)];
        else
            stateR(:,s) = [swg_traj(ini(n+1)-1+q,2); swg_traj(ini(n+1)-1+q,1); swg_traj(ini(n+1)-1+q,3)];
            stateL(:,s) = [u0y; u0x; -0];
        end
    end
%     animate(stateR(:,s), stateL(:,s), stateCoM(:,s), animateOn, robot, hLeftRel, hRightRel, hCoMRel, speedupfactor, s, z_robot);    

end
%% plot result
% r = 0.07;
% hold on
% [x,y,z] = sphere(50);
% x0 = stateCoM(1,end); y0 = stateCoM(2,end); z0 = z_robot+0.2;
% x = x*r + x0;
% y = y*r + y0;
% z = z*r + z0;
% surface(x,y,z,'FaceColor', 'b','EdgeColor','none')
% plot3([x0 x0],[y0 y0],[z0 z_robot],'LineWidth',1,'Color','b');
% %%
% r = 0.05;
% hold on
% [x,y,z] = sphere(50);
% x0 = stateL(1,end); y0 = stateL(2,end); z0 = stateL(3,end);
% x = x*r + x0;
% y = y*r + y0;
% z = z*r + z0;
% surface(x,y,z,'FaceColor', 'b','EdgeColor','none')
% %%
% r = 0.05;
% hold on
% [x,y,z] = sphere(50);
% x0 = stateR(1,end); y0 = stateR(2,end); z0 = stateR(3,end);
% x = x*r + x0;
% y = y*r + y0;
% z = z*r + z0;
% surface(x,y,z,'FaceColor', 'b','EdgeColor','none')
% %%
% xlabel('Y(m)');
% ylabel('X(m)');
% zlabel('Z(m)');
% hold off
figure(2)
plot(XI_ref_X(1,:),XI_ref_X(2,:),'color','k','LineStyle','--','linewidth',1.5);hold on;
plot(ZETA_mea_x_3Mass(1,:),ZETA_mea_x_3Mass(2,:),'color','k','linewidth',2);hold on;
plot(CoMx_3Mass(1,:),x_com(2:end,1),'color','g','LineStyle','--','linewidth',1.5);hold on;
plot(CoMx_3Mass(1,:),CoMx_3Mass(2,:),'color','g','linewidth',2);hold on;
plot(U0_x(1,:),U0_x(2,:),'color','c','linewidth',1.5);hold on;
plot(UT_x(1,:),UT_x(2,:),'color','b','linewidth',1.5);hold on;
plot(ZETA_mea_x_3Mass(1,:),PcZMP_X,'color','r','linewidth',2);
plot(ZMP_PEND(1,:), ZMP_PEND(2,:),'color','m','LineStyle','-','linewidth',1.5)
legend('\xi_{ref,x}','\xi_{meas,x}','x_{com,ref}','x_{com,meas}','u_{0,x}','u_{T,x}','P_{cZMP,x} + u_{0,x}','zmp_{pend}') 
xlabel('time(s)');
ylabel('position_{x} (m)');
grid on

%%
figure(3)
plot(XI_ref_Y(1,:),XI_ref_Y(2,:),'color','k','LineStyle','--','linewidth',1.5);hold on;
plot(ZETA_mea_y_3Mass(1,:),ZETA_mea_y_3Mass(2,:),'color','k','linewidth',2);hold on;
plot(CoMy_3Mass(1,:),x_com(2:end,2),'color','g','LineStyle','--','linewidth',1.5);hold on;
plot(CoMy_3Mass(1,:),CoMy_3Mass(2,:),'color','g','linewidth',2);hold on;
plot(U0_y(1,:),U0_y(2,:),'color','c','linewidth',1.5);
plot(UT_y(1,:),UT_y(2,:),'color','b','linewidth',1.5);hold on;
plot(ZETA_mea_y_3Mass(1,:),PcZMP_Y,'color','r','linewidth',2);
plot(ZMP_PEND(1,:), ZMP_PEND(3,:),'color','m','LineStyle','-','linewidth',1.5)
legend('\xi_{ref,y}','\xi_{meas,y}','y_{com,ref}','y_{com,meas}','u_{0,y}','u_{T,y}','P_{cZMP,y} + u_{0,y}','zmp_{pend}')
xlabel('time(s)');
ylabel('position_{y} (m)');
grid on

% figure(4)
% plot(ZETA_err_x(1,:),ZETA_err_x(2,:),'color','k','LineStyle','--','linewidth',2);hold on;
% plot(ZETA_mea_x_3Mass(1,:),PcZMP_XX,'color','r','linewidth',2);
% plot(UT_xEr(1,:),UT_xEr(2,:),'color','b','linewidth',2);hold on;
% plot(bT_xEr(1,:),bT_xEr(2,:),'color','m','linewidth',2);hold on;
% legend('\xi_{err,x}','P_{cZMP,x}','u_{T,err}','b_{T,err}')
% xlabel('time(s)');
% ylabel('distance_{x} (m)');
% %ylim([-0.09 0.09])
% grid on
% 
% figure(5)
% plot(ZETA_err_y(1,:),ZETA_err_y(2,:),'color','k','LineStyle','--','linewidth',2);hold on;
% plot(ZETA_mea_y_3Mass(1,:),PcZMP_YY,'color','r','linewidth',2);
% plot(UT_yEr(1,:),UT_yEr(2,:),'color','b','linewidth',2);hold on;
% plot(bT_yEr(1,:),bT_yEr(2,:),'color','m','linewidth',2);hold on;
% legend('\xi_{err,y}','P_{cZMP,y}','u_{T,err}','b_{T,err}')
% xlabel('time(s)');
% ylabel('distance_{y} (m)');
% %ylim([-0.05 0.05])
% grid on
% 
% figure(6)
% plot(T_step(1,:),T_step(2,:),'color','r','linewidth',2);hold on;
% plot(t_var(1,:),t_var(2,:),'color','k','LineStyle','--','linewidth',2);hold on;
% legend('T_{new}','t_{curr}')
% xlabel('time(s)');
% ylabel('time(s)');
% xlim([-0.1 6])
% set(gca, 'DataAspectRatio',[5 1 1])
% grid on

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
function animate(stateR, stateL, stateCoM, animateOn, robot, hLeftRel, hRightRel, hCoMRel, speedupfactor, idx, z_robot)

n = [0;  0; -1]; % x
s = [-1; 0; 0];  % y
a = [0;  1; 0];  % z
R = [n s a];   

% Get Left joints
p = stateR-stateCoM.*[1; 1; 1];
transmatL =  [R     p; 
            [0 0 0 1]];
isLeft = true; 
qLeft = invKinBody2Foot(transmatL, isLeft); % Call IK function

% Get Right joints
p = stateL-stateCoM.*[1; 1; 1]; 
transmatR =  [R     p; 
            [0 0 0 1]];
isLeft = false; 
qRight = invKinBody2Foot(transmatR, isLeft);
appendLine(hLeftRel, stateL); 
appendLine(hRightRel, stateR);
appendLine(hCoMRel, stateCoM-[0;0;-0.2]);
% updateLine(hLeftPoint, [stepinfos{footidx}.state([1,3,5],sidx) stateL([1,3,5],sidx)]); 
% updateLine(hRightPoint, [stepinfos{footidx}.state([1,3,5],sidx) stateR([1,3,5],sidx)]); 


% Animate
if animateOn
    if rem(idx,speedupfactor) == 0
        updateJoints(robot, qRight, qLeft, stateCoM);
    end 
end  
end
function [robot, hLeftRel, hRightRel, hCoMRel] = createRobot(stateC, z_robot, stateL, stateR, stateCoM)

% NOTE: make sure parameters match in inverse kinematics function 
L1 = 0.12; 
L2 = 0; 
L3 = 0.4;
L4 = 0.38;
L5 = 0;
       
robot = rigidBodyTree;
% Right Leg
dhparams = [0       0        0      0;
            0       0        0      0;     % Base -> pelvisy
            0       0        0      0;     % pelvisy -> pelvisx
            L1      0       -L2     0;     % pelvisx -> hip yaw
            0      -pi/2     0      0;     % Hip yaw -> hip roll
            0      +pi/2     0      0;     % Hip roll -> hip pitch       
            L3      0        0      0;     % Hip pitch -> knee pitch
            L4      0        0      0;     % Knee pitch -> ankle pitch
            0       pi/2     0      0;     % Ankle pitch -> ankle roll
            L5      0        0      0];    % Ankle roll -> end effector (foot)
       
for idx = 1:size(dhparams,1)
    rightLeg(idx) = rigidBody("rightleg"+idx);
    rightJnt(idx) = rigidBodyJoint("rightjnt"+idx, 'revolute');
    if idx==1 || idx==2 || idx==3
        rightJnt(idx) = rigidBodyJoint("rightjnt"+idx, 'prismatic');
        if idx==1
            rightJnt(idx).JointAxis = [0 0 1];
        elseif idx==2
            rightJnt(idx).JointAxis = [0 -1 0];
        else
           rightJnt(idx).JointAxis = [1 0 0];
        end
    end
    setFixedTransform(rightJnt(idx),dhparams(idx,:),'dh'); 
    rightLeg(idx).Joint = rightJnt(idx);
    if idx==1
        addBody(robot,rightLeg(idx),"base");     
    else
        addBody(robot,rightLeg(idx),"rightleg"+(idx-1));
    end
    
end

% Left Leg
dhparams = [0      0         0      0;
            0      0         0      0;
            0      0         0      0;
           -L1     0        -L2     0;    % Only difference with right leg is the
            0     -pi/2      0      0;    % first element is -L1 instead of L1
            0     +pi/2      0      0; 
            L3      0        0      0;
            L4      0        0      0;
            0       pi/2     0      0;
            L5      0        0      0];
        
for idx = 1:size(dhparams,1)
    leftLeg(idx) = rigidBody("leftleg"+idx); 
    leftJnt(idx) = rigidBodyJoint("leftjnt"+idx, 'revolute');
    if idx==1 || idx==2 || idx==3
        leftJnt(idx) = rigidBodyJoint("leftjnt"+idx, 'prismatic');
        if idx==1
            leftJnt(idx).JointAxis = [0 0 1];
        elseif idx==2
            leftJnt(idx).JointAxis = [0 -1 0];
        else
            leftJnt(idx).JointAxis = [1 0 0];
        end
    end
    setFixedTransform(leftJnt(idx),dhparams(idx,:),'dh'); 
    leftLeg(idx).Joint = leftJnt(idx);
    if idx==1
        addBody(robot,leftLeg(idx),"base");     
    else
        addBody(robot,leftLeg(idx),"leftleg"+(idx-1));
    end
end

% showdetails(robot)
hFig = figure; 
% hFig.Visible = 'on'; 
hFig.Units = 'Normalized'; 
hFig.OuterPosition = [0 0 1 1];
% hFig.Position = [0.13 0.13 400 400];
hAx2 = axes(hFig);
mcolors = get(gca, 'colororder'); 

desconfig = robot.homeConfiguration;


qright0 = zeros(1,6); 
qleft0 = zeros(1,6);
updateJoints(robot, qright0, qleft0, stateC)
hold on
hLeftRel = plot3(stateL(1),stateL(2),stateL(3),'Color',mcolors(1,:),"LineWidth", 2); 
hRightRel = plot3(stateR(1),stateR(2),stateR(3),'Color',mcolors(2,:),"LineWidth", 2);
hCoMRel = plot3(stateCoM(1),stateCoM(2),stateCoM(3) + 0.2,'Color',mcolors(3,:),"LineWidth", 2.5);

view(3)
grid off
axis([-1.1*z_robot 0.5*z_robot -0.5*z_robot 3.5*z_robot  -.1*z_robot 1.5*z_robot])

end
function updateJoints(robot, anglesright, anglesleft, stateC)
    desconfig = robot.homeConfiguration;
    
    desconfig(1).JointPosition = 0.68;
    desconfig(2).JointPosition = -stateC(2); % angle offset
    desconfig(3).JointPosition = stateC(1);
    for idx = 1:length(anglesright)
        desconfig(idx+4).JointPosition = anglesright(idx);
    end 
    desconfig(5).JointPosition = desconfig(5).JointPosition - pi;
    desconfig(6).JointPosition = desconfig(6).JointPosition + pi/2; 
    
    desconfig(11).JointPosition = 0.68;
    desconfig(12).JointPosition = -stateC(2);
    desconfig(13).JointPosition = stateC(1);
    for idx = 1:length(anglesleft)
        desconfig(idx+14).JointPosition = anglesleft(idx);
    end 
    desconfig(15).JointPosition = desconfig(15).JointPosition - pi; 
    desconfig(16).JointPosition = desconfig(16).JointPosition + pi/2; 
       
    % update graphics 
    show(robot, desconfig, 'PreservePlot', false);
%     title('Walking Pattern Inverse Kinematics')
    pause(0.001)
end
function appendLine(gHandle, points)
    gHandle.XData(end+1) = points(1); 
    gHandle.YData(end+1) = points(2); 
    gHandle.ZData(end+1) = points(3); 
end
function updateLine(gHandle, points)
    gHandle.XData = points(1,:); 
    gHandle.YData = points(2,:); 
    gHandle.ZData = points(3,:); 
end