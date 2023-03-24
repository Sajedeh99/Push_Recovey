function [zmp_pend, xi_ini, xi_eos] = input3Mass(is_left, Lp, Wnom, N, Lnom, delta_z_vrp, swingHeight, T, t_sample, mswg, msup, mpend, mfeet, m)
g = 9.807;
omega = sqrt(g/delta_z_vrp);
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
for ith = 1:N+2
    indx =  (T*(ith-1)/t_sample);
        
    endtime = T/t_sample;
    init_time = indx*t_sample + t_sample;
    final_time = init_time + T - t_sample;
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
    swing_traj(indx+1:indx+endtime,:) = [qswing', dqswing', ddqswing'];
    for j = 1:endtime
        p_ref(int32(indx+j),:) = r_vrp(ith,:);
        
        % X direction
        Mfeet(int32(indx+j),1) = mswg*(qswing(1,j) - p_ref(int32(indx+j),1))*(g + ddqswing(3,j)) ...
            - mswg*(ddqswing(1,j))*(qswing(3,j) - delta_z_vrp);
        Mfeet(int32(indx+j),2) = mswg*(qswing(2,j) - p_ref(int32(indx+j),2))*(g + ddqswing(3,j)) ...
            - mswg*(ddqswing(2,j))*(qswing(3,j) - delta_z_vrp);
        
        zmp_feet(int32(indx+j),1) = Mfeet(int32(indx+j),1)/(mfeet*g) + p_ref(int32(indx+j),1);
        zmp_feet(int32(indx+j),2) = Mfeet(int32(indx+j),2)/(mfeet*g) + p_ref(int32(indx+j),2);
        
        zmp_pend(int32(indx+j),1) = (m/mpend)*p_ref(int32(indx+j),1) - (mfeet/mpend)*zmp_feet(int32(indx+j),1);
        zmp_pend(int32(indx+j),2) = (m/mpend)*p_ref(int32(indx+j),2) - (mfeet/mpend)*zmp_feet(int32(indx+j),2);
    end
end
zmp_pend(:,3) = delta_z_vrp*ones(length(zmp_pend),1);

[xi_ini, xi_eos] = Xi(N, zmp_pend, omega, T, t_sample);

function [xi_ini, xi_eos] = Xi(N, r_vrp, omega, t_step, t_sample)
xi_eos = zeros(N+2,3);
xi_eos(N+2,:) = r_vrp(end,:);
for k = N+1:-1:1
        xi_eos(k,:) = r_vrp((k+1)*int32(t_step/t_sample),:) + (exp(-omega*t_step))*(xi_eos(k+1,:)-r_vrp((k+1)*int32(t_step/t_sample),:));
        xi_ini(k+1,:) = xi_eos(k,:);
end
xi_ini(1,:) = r_vrp(k*int32(t_step/t_sample),:) + (exp(-omega*t_step))*(xi_eos(k,:)-r_vrp(k*int32(t_step/t_sample),:));
end
end