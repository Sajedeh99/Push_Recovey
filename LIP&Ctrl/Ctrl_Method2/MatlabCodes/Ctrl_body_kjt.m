clear all; clc; close all;
N = 3;
is_left = false;

global t_sample
t_sample = 0.001;

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

t_ini = 2;
t_fnl = 2;
NL = 1.6/t_sample;

[Tnom,Lnom,Wnom,tau_nom] = Param_fcn2(L_min0, L_max0, W_min0, W_max0, T_min, T_max, omega, Vx, Vy);


traj_ds = struct('index', {}, 'xi_ds', {}, 'state', {});
traj_ss = struct('index', {}, 'xi_ss', {}, 'state', {});

for i=1:N+3
    if i == 1
        r_f_r(i,1) = 0;
        r_f_l(i,1) = 0;
    elseif i == N+3
        r_f_l(i, 1) = (Lnom)*(i-2);
        r_f_r(i,1) = (Lnom)*(i-2);
    elseif mod(i,2) == 0
        r_f_l(i,1) = (Lnom)*(i-1);
        r_f_r(i,1) = r_f_r(i-1,1);  
    else
        r_f_r(i,1) = (Lnom)*(i-1);
        r_f_l(i,1) = r_f_l(i-1,1);
    end
    r_f_l(i,2:3) = [(Lp/2 + Wnom) 0];
    r_f_r(i,2:3) = [-(Lp/2 + Wnom) 0];
end
ini_org(1) = 1; fnl_org(1) = Tnom/t_sample;
for i = 2:N+3
    ini_org(i) = (i-1)*int32(Tnom/t_sample) + (i-1);
    fnl_org(i) = ini_org(i) + int32(Tnom/t_sample);
    if i==N+3
        fnl_org(i) = ini_org(i) + int32(t_fnl/t_sample) - 1;
    end
end
ini_org(N+4) = fnl_org(N+3)+1;
r_vrp = r_f_r;
r_vrp(2:2:N+3, 1:2) = r_f_l(2:2:N+3,1:2);
if is_left == false
    temp = r_f_l;
    r_f_l = [1 -1 1].*r_f_r;
    r_f_r = [1 -1 1].*temp;
    r_vrp = r_f_l;
    r_vrp(2:2:N+3, 1:2) = r_f_r(2:2:N+3,1:2);
end
r_vrp(:,3) =+ delta_z_vrp;

[xi_ini, xi_eos] = Xi(N, r_vrp, omega, Tnom);
for ith = 1:N+2
    b_nom(ith,:) = (xi_eos(ith,:) - r_vrp(ith+1,:));
end

incrs = int32(Tnom/t_sample);
p_ref = zeros(t_ini/t_sample,3);
for n_ = 1:3
    for t_ = 1:int32(Tnom/t_sample)
        p_ref(int32(t_ini/t_sample)+ (n_-1)*incrs + t_,:) = r_vrp(n_,:);
    end
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
Step = 1; i = 1; q = 0; n = 0; s = 1; m = 0;
qpresult = [0;0;0;0;0]; inc(s)=0;
ini = zeros(1,N+4); ini(1)=1; fnl = zeros(1,N+3);

%% preview Ctrl parameters
A = [1 t_sample (t_sample^2)/2;
     0    1       t_sample;
     0    0           1];
B = [(t_sample^3)/6; (t_sample^2)/2; t_sample];
C = [1 0 -delta_z_vrp/g];

Qe = 1; R = 1e-6;
Qx = zeros(3,3);
B_ = [C*B; B];
I_ = [1; zeros(3,1)];
F_ = [C*A; A];
Q_ = [Qe zeros(1,3); zeros(3,1) Qx];
A_ = [I_ F_];

[K_,L,G,info] = idare(A_,B_,Q_,R,[],[]);

GI = (B_'*K_*B_ + R)\(B_'*K_*I_);
Gx = (B_'*K_*B_ + R)\(B_'*K_*F_);
Ac_ = A_-B_*inv(R+B_'*K_*B_)*B_'*K_*A_;

for l = 1:NL
    if l == 1;
        X_(:,l) = -Ac_'*K_*I_;
    else
        X_(:,l) = Ac_'*X_(:,l-1);
    end
end

CoM = struct('x', zeros(3,1), 'y', zeros(3,1));
p_out = zeros(1,2);
u_ref = zeros(1,2);
indx = 0;

%% COM generation for t_ini
    for k = 1:t_ini/t_sample             
        % ref
        err_ref = [0 0];
        p_out(k+indx,:) = C*[CoM.x(:,k+indx) CoM.y(:,k+indx)];
        
        DCM.x(k+indx) = CoM.x(1,k+indx) + CoM.x(2,k+indx)/omega;
        DCM.y(k+indx) = CoM.y(1,k+indx) + CoM.y(2,k+indx)/omega;
        
        for i_ = 1:k+indx
            err_ref = err_ref + p_out(i_,:)-(p_ref(i_,1:2));
        end
        preview_term = [0 0];
        for j = 1:NL
            if j==1
                Gp(j) = -GI;
            else
                Gp(j) = (B_'*K_*B_ + R)\(B_'*X_(:,j-1));
            end
            preview_term = preview_term + Gp(j)*p_ref(k+indx+j,1:2);
        end
        u_ref(k+indx,:) = -GI*err_ref - Gx*[CoM.x(:,k+indx) CoM.y(:,k+indx)] - preview_term;
        CoM.x(:,k+indx+1) = A*CoM.x(:,k+indx) + B*u_ref(k+indx,1);
        CoM.y(:,k+indx+1) = A*CoM.y(:,k+indx) + B*u_ref(k+indx,2);
    end
indx = int32(t_ini/t_sample); 
    
%% control loop
while Step(i) == 1
    q = q+1;
    s = s + 1;
    m = m + 1;
    % update pattern parameter
    for j = n:N+1
       r_vrp(j+2,1) = r_vrp(j+2,1) + qpresult(1);
       r_vrp(j+2,2) = r_vrp(j+2,2) + qpresult(3);
    end
    
    % p_ref generation
    for ith = n+1:N+3
        inc(s) = floor(T/t_sample) - int32(Tnom/t_sample);
        ini(ith+1:end) = ini_org(ith+1:end) + inc(s).*ones(1,N+4-ith);
        fnl(ith:end) = fnl_org(ith:end) + inc(s).*ones(1,N+4-ith);
        p_ref(int32(t_ini/t_sample)+ini(ith):int32(t_ini/t_sample)+fnl(ith),:) = r_vrp(ith,:).*ones(fnl(ith)-ini(ith)+1,1);
    end
    
    %% CoM generation
    err_ref = [0 0];
    p_out(m+indx,:) = C*[CoM.x(:,m+indx) CoM.y(:,m+indx)];
    
    DCM.x(m+indx) = CoM.x(1,m+indx) + CoM.x(2,m+indx)/omega;
    DCM.y(m+indx) = CoM.y(1,m+indx) + CoM.y(2,m+indx)/omega;
    for i_ = 1:m+indx
        err_ref = err_ref + p_out(i_,:)-(p_ref(i_,1:2));
    end    
    preview_term = [0 0];
    for j = 1:NL
        if j==1
            Gp(j) = -GI;
        else
            Gp(j) = (B_'*K_*B_ + R)\(B_'*X_(:,j-1));
        end
        preview_term = preview_term + Gp(j)*p_ref(m+indx+j,1:2);
    end
    u_ref(m+indx,:) = -GI*err_ref - Gx*[CoM.x(:,m+indx) CoM.y(:,m+indx)] - preview_term;
    CoM.x(:,m+indx+1) = A*CoM.x(:,m+indx) + B*u_ref(m+indx,1);
    CoM.y(:,m+indx+1) = A*CoM.y(:,m+indx) + B*u_ref(m+indx,2);
   %%
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
    [qpresult, Opt_Vector] = controller2(t, T, Lnom, Wnom, L_min, L_max, W_min, W_max, T_min, T_max,...
          b_nom(n+1,1), b_nom(n+1,2), omega, zeta_mea_x, zeta_mea_y, r_vrp(n+2,1), r_vrp(n+2,2),...
          zeta_err_x, zeta_err_y, PcZMP_y(q,n+1), PcZMP_x(q,n+1)); %PcZMP_y(q,n+1), PcZMP_x(q,n+1)

    T = (1/omega)*log(Opt_Vector(3));

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
    
    [Opt_Vector(1); Opt_Vector(2); Opt_Vector(3); Opt_Vector(4); Opt_Vector(5)]
    [n T t]
    
    % going next step
    if t>=T
        t = 0;
        Ts(i) = T;
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
        ini_org = ini;
        fnl_org = fnl;
       
    end

    if n == N+2 % N
        break
    end
end

%% plot result
t_sim = t_sample:t_sample:length(p_ref)*t_sample;
figure(1)
plot(ZETA_mea_x(1,:),ZETA_mea_x(2,:),'color','g');hold on;
plot(XI_ref_X(1,:),XI_ref_X(2,:),'color','k','LineStyle','-','linewidth',2);hold on;
plot(CoMx(1,:),CoMx(2,:),'color','m');hold on;
plot(UT_x(1,:),UT_x(2,:),'color','b');hold on;
plot(U0_x(1,:),U0_x(2,:),'color','c','linewidth',2);hold on;
figure(3)
plot(t_sim,p_ref(:,1),'color','r','LineStyle',':','linewidth',2);hold on;
plot(t_sim(1:end-int32(t_fnl/t_sample)+1),CoM.x(1,:),'color','m','LineStyle',':','linewidth',2);hold on;
plot(t_sim(1:end-int32(t_fnl/t_sample)),DCM.x(1,:),'color','k','LineStyle',':','linewidth',2);hold on;
% plot(ZETA_mea_x(1,:),PcZMP_X,'color','r','linewidth',2);
figure(2)
plot(ZETA_mea_y(1,:),ZETA_mea_y(2,:),'color','g','linewidth',2);hold on;
plot(XI_ref_Y(1,:),XI_ref_Y(2,:),'color','k','LineStyle','-','linewidth',2);hold on;
plot(CoMy(1,:),CoMy(2,:),'color','m','linewidth',2);hold on;
plot(UT_y(1,:),UT_y(2,:),'color','b','linewidth',2);hold on;
plot(U0_y(1,:),U0_y(2,:),'color','c','linewidth',2);hold on;
figure(4)
plot(t_sim,p_ref(:,2),'color','r','LineStyle',':','linewidth',2);hold on;
plot(t_sim(1:end-int32(t_fnl/t_sample)+1),CoM.y(1,:),'color','m','LineStyle',':','linewidth',2);hold on;
plot(t_sim(1:end-int32(t_fnl/t_sample)),DCM.y(1,:),'color','k','LineStyle',':','linewidth',2);hold on;
% plot(ZETA_mea_y(1,:),PcZMP_Y,'color','r','linewidth',2);

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
