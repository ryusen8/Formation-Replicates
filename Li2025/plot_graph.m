G = graph(adj_mat);
adj_mat_new = omega_new - diag(diag(omega_new));
G_new = graph(adj_mat_new);
figure;
subplot(1,2,1);
p = plot(G);
p.XData = nominal_config(:,1);
p.YData = nominal_config(:,2);
p.LineWidth = 1.5;
p.MarkerSize = 8;
subplot(1,2,2);
p_new = plot(G_new); hold on
p_new.XData = nominal_config_new(:,1);
p_new.YData = nominal_config_new(:,2);
p_new.LineWidth = 1.5;
p_new.MarkerSize = 8;