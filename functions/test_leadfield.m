function test_leadfield(leadfield, brain_mesh,elec,source_index)
% function to test the leadfield
% from https://www.fieldtriptoolbox.org/workshop/oslo2019/forward_modeling/

figure('units', 'normalized', 'outerposition', [0 0 0.5 0.5])
% source_index: a superficial sources
sensory_dipole_current = 100e-9; % Am (realistic)

n_sensors = length(elec.label);

inside_sources = find(leadfield.inside);
inside_index = inside_sources(source_index);
lead = leadfield.leadfield{inside_index};
xs = zeros(1, n_sensors);
ys = zeros(1, n_sensors);
zs = zeros(1, n_sensors);
voltages = zeros(1, n_sensors);
titles = {'Lead field (x)' 'Lead field (y)' 'Lead field (z)'};

% get the xyz and norm

for sensor_index = 1:n_sensors
    this_x = lead(sensor_index, 1);
    this_y = lead(sensor_index, 2);
    this_z = lead(sensor_index, 3);
    this_norm = norm(lead(sensor_index, :));
    xs(sensor_index) = this_x * sensory_dipole_current;
    ys(sensor_index) = this_y * sensory_dipole_current;
    zs(sensor_index) = this_z * sensory_dipole_current;
    voltages(sensor_index) = this_norm * sensory_dipole_current;
end

% plot xyz
axes = {xs ys zs};

for axis_index = 1:3
    this_axis = axes{axis_index};
    subplot(1, 3, axis_index)
    hold on
    ft_plot_topo3d(elec.chanpos, this_axis, 'facealpha', 0.7)
    c = colorbar('location', 'southoutside');
    c.Label.String = 'Lead field (V)';
    axis tight
    ft_plot_mesh(brain_mesh, 'facealpha', 0.10);
    ft_plot_sens(elec, 'elecsize', 20);
    title(titles{axis_index})
    src_pos = [leadfield.pos(inside_index, 1), ...
        leadfield.pos(inside_index, 2), ...
        leadfield.pos(inside_index, 3)];
    src_ori = 1*circshift([1 0 0], axis_index-1);
%     plot3(leadfield.pos(inside_index, 1), ...
%         leadfield.pos(inside_index, 2), ...
%         leadfield.pos(inside_index, 3), 'bo', ...
%         'markersize', 20, 'markerfacecolor', 'r')
    quiver3(src_pos(1), src_pos(2), src_pos(3), src_ori(1), src_ori(2),...
        src_ori(3), 25, 'LineWidth', 5, 'Color', 'r', 'MaxHeadSize', 30)
    caxis(max(abs(this_axis))*[-1 1])
end

%% plot norm
% 
% figure('units', 'normalized', 'outerposition', [0 0 0.5 0.85])
% hold on
% ft_plot_topo3d(elec.chanpos, voltages, 'facealpha', 0.8)
% if strcmp(headmodel.type, 'dipoli')
%     caxis([0, 10e-6])
% end
% c = colorbar('location', 'eastoutside');
% c.Label.String = 'Lead field (V)';
% axis tight
% ft_plot_mesh(headmodel.bnd(3), 'facealpha', 0.10);
% ft_plot_sens(elec, 'elecsize', 20);
% title('Leadfield magnitude')
% plot3(leadfield.pos(inside_index, 1), ...
%   leadfield.pos(inside_index, 2), ...
%   leadfield.pos(inside_index, 3), 'bo', ...
%   'markersize', 20, 'markerfacecolor', 'r')
% 
% view(-90, 0)