function omega_add_node = omega_add_node(agents_num, perceived_node, nominal_config, add_config, follower_pos, add_pos)
%{
perceived_node: 3 x 1 | [i j k] | i<j<k
choose_node: 4 x 1 | [i j k u] | i<j<k<u
nominal_config: agents_num x 2
nominal_config_new: (agents_num+1) x 2
follower_pos: follower_num x 2
follower_pos_new: (follower_num+1) x 2
omega_add_vertex: (agents_num+1) x (agents_num+1)
%}
    scaling_param = 1;
    choose_node = [perceived_node, agents_num+1];
    agents_index = linspace(1,agents_num,agents_num); % 未加点前的agent编号
    new_order = [setdiff(agents_idx, choose_node), choose_node]; % 重新排列编号
    nominal_config_new = [nominal_config; add_config];
    follower_pos_new = vertcat(follower_pos, add_pos);
    % 解方程Pu*φ = 0
    Pu_1 = nominal_config_new(choose_node,1)';
    Pu_2 = nominal_config_new(choose_node,2)';
    Pu_3 = ones(1,4);
    Pu = vertcat(Pu_1,Pu_2,Pu_3);
    phi = null(Pu,'rational');
    omega_u = scaling_param * (phi * phi');
    omega_a = [
        omega,             zeros(agent_nums,1);
        zeros(1,agent_nums),     0                ];
    omega_b = [
        zeros(agent_nums-3),         zeros(agent_nums-3, 4);
        zeros(4,agent_nums-3),       omega_u          ];
    omega_add_node = omega_a + omega_b;

