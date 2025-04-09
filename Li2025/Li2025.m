% LI H, CHEN H, WANG X, 等. Distributed Framework Construction for Affine Formation Control[EB/OL]. 
% arXiv, 2025[2025-03-10]. http://arxiv.org/abs/2501.01817.
clear;
clc;
% 仿真步长
dt = 0.01;
sim_time = 40;
sim_steps = sim_time/dt;

% 维数和个数
dim = 2;
agent_num = 7;
leader_num = 3;
follower_num = agent_num - leader_num;

% 图矩阵
adj_mat = [
    0       0.2741      0.2741      -0.137      -0.137      0      0     ;
    0.2741  0           0            0.5482      0          0     -0.137 ;
    0.2741  0           0            0           0.5482    -0.137  0     ;
    -0.137  0.5482      0            0           0.0685     0.2741 0     ;
    -0.137  0           0.5482       0.0685      0          0      0.2741;
    0       0          -0.137        0.2741      0          0      0.137 ;
    0       -0.137      0            0           0.2741     0.137  0     ;];
degree_mat = diag(sum(adj_mat,1));
stress_mat = degree_mat - adj_mat;

% 坐标
agent_traj_x = zeros(sim_steps, agent_num+1);
agent_traj_y = zeros(sim_steps, agent_num+1);
nominal_config = [2  0; 1  1; 1 -1; 0  1; 0 -1; -1  1; -1 -1];
leader_pos = nominal_config(1:leader_num, :);
follower_pos_0 = randi(10,follower_num,dim);
follower_pos = follower_pos_0;
Kp = 0.5;
Ki = 1;
% leader速度
leader_vel = [1, 0];
% 跟踪误差
tracking_error = zeros(sim_steps, follower_num);

% t = 20s 时加入新点
add_node_time = 20;
scaling_param = 1;

% 加点前控制编队达到应力平衡
step = 1;
error_0 = zeros(dim*follower_num, 1);
% pos_error = zeros(follower_num, dim);
while(step*dt < add_node_time)
    leader_pos = leader_pos + leader_vel * dt;
    [follower_vel, error] = follower_input(Kp, Ki, dim, dt, stress_mat, leader_num, leader_pos, error_0);
    error_0 = error;
    follower_pos = follower_pos + follower_vel * dt;
    agent_x = vertcat(leader_pos(:,1),follower_pos(:,1))';
    agent_y = vertcat(leader_pos(:,2),follower_pos(:,2))';
    agent_traj_x(step,1:agent_num) = agent_x;
    agent_traj_y(step,1:agent_num) = agent_y;
    desired_follower_pos = nominal_config(leader_num+1:end,:) + leader_vel;
    tracking_error(step, :) = vecnorm(follower_pos - desired_follower_pos, 2, 2);
    step = step + 1;
end

% 加入一个新点 node 8
follower_num = follower_num + 1;
add_pos = [0 2];
follower_pos = vertcat(follower_pos, add_pos);
nominal_config_new = [nominal_config; -2, 0];
tracking_error = [tracking_error, zeros(sim_steps, 1)];

% 计算加点后的应力矩阵 (这里应该用nominal_config来计算而不是follower_pos)
% 选择三个点后重新排列index
% 原应力矩阵要按排列后的index改写，这里选的是末尾三个所以应力矩阵和原来一样
choose_node = [5,6,7,8]; % i<j<k<u
agent_idx = linspace(1,agent_num,agent_num);
new_order = [setdiff(agent_idx, choose_node), choose_node];

% 解方程Puφ = 0
Pu_1 = nominal_config_new(choose_node,1)';
Pu_2 = nominal_config_new(choose_node,2)';
Pu_3 = ones(1,4);
Pu = vertcat(Pu_1,Pu_2,Pu_3);
phi = null(Pu,'rational');

omega_u = scaling_param * (phi * phi');
omega_a = [
    stress_mat,             zeros(agent_num,1);
    zeros(1,agent_num),     0                ];
omega_b = [
    zeros(agent_num-3),         zeros(agent_num-3, 4);
    zeros(4,agent_num-3),       omega_u          ];
omega_new = omega_a + omega_b;
isUniversallyRigid(omega_new, 2);

agent_num = agent_num + 1;
error_0 = vertcat(error_0, zeros(2,1));
% disp(step);
% 用新应力矩阵控制
while(step*dt >= add_node_time && step*dt <= sim_time)
    leader_pos = leader_pos + leader_vel * dt;
    [follower_vel, error] = follower_input(Kp, Ki, dim, dt, omega_new, leader_num, leader_pos, error_0);
    error_0 = error;
    follower_pos = follower_pos + follower_vel * dt;
    agent_x = vertcat(leader_pos(:,1),follower_pos(:,1))';
    agent_y = vertcat(leader_pos(:,2),follower_pos(:,2))';
    agent_traj_x(step,1:agent_num) = agent_x;
    agent_traj_y(step,1:agent_num) = agent_y;
    desired_follower_pos = nominal_config_new(leader_num+1:end,:) + leader_vel;
    tracking_error(step, :) = vecnorm(follower_pos - desired_follower_pos, 2, 2);
    step = step + 1;    
end

plot(agent_traj_x, agent_traj_y); hold on; 

% 先不考虑仿射变换，只平移
% 维数有点问题
function [follower_vel, error] = follower_input(Kp, Ki, dim, dt, stress_mat, leader_num, leader_pos, error_0)
    % leader_pos nl*d   |   follower_pos nf*d
    omega_bar = kron(stress_mat, eye(dim)); % dn*dn
    omega_bar_ff = omega_bar(dim*(leader_num)+1:end, dim*(leader_num)+1:end); % dnf*dnf
    omega_bar_fl = omega_bar(dim*(leader_num)+1:end, 1:dim*leader_num); % dnf*dnl
    p_l = pos2vec(leader_pos); % dnl
    p_f = -inv(omega_bar_ff) * omega_bar_fl * p_l; % dnf
    d_error = omega_bar_ff * p_f + omega_bar_fl * p_l; % dnf
    error = error_0 + d_error * dt; % dnf
    u_f = -Kp * omega_bar_ff * p_f - Kp * omega_bar_fl * p_l - Ki * error;
    % follower_vel = u_f(1:dim:end, 1:dim:end);
    follower_vel = vec2pos(u_f);
end

function vec = pos2vec(pos)
    vec = reshape(pos', numel(pos), 1);
end

function pos = vec2pos(vec)
    pos = reshape(vec, 2, length(vec)/2)';
end