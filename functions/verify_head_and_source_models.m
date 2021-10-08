function fig = verify_head_and_source_models(mesh, sourcemodel, elec, ortho_view)

% draw 3D representation of headmodel, source points and electrodes
fig = figure('Visible', 'off');
ft_plot_mesh(mesh(1), 'facecolor',[0.9 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(mesh(2),'edgecolor','none','facealpha',0.4);
ft_plot_mesh(mesh(3),'edgecolor','none','facecolor',0.4*[1 1 1],'facealpha', 0.3);
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:),'vertexcolor', 'k');
if exist('elec', 'var') && ~isempty(elec)
    ft_plot_sens(elec, 'style', '.b', 'label', 'label')
end
hold off

if ~exist('ortho_view', 'var') || ~ortho_view
    return
end

view_points = [0 0; 90 0; 0 90];
% view_points = [-180 0; 90 0; -180 90];
figure('color', 'w', 'Position', [2, 42, 798, 774],'Visible', 'on');
ax(1) = subplot(2,2,1);
ax(2) = subplot(2,2,2);
ax(3) = subplot(2,2,3);
arrayfun(@(axx) copyobj(fig.Children.Children,axx), ax)
arrayfun(@(idx) view(ax(idx), view_points(idx,:)), 1:3)
axis(ax, 'vis3d', 'equal', 'off')
ax(1).Position = [0,0.5,0.48,0.48];
ax(2).Position = [0.5,0.5,0.48,0.48];
ax(3).Position = [0,0.0,0.48,0.48];
text(ax(3),0.1,0.90,'L','Units','normalized', 'Fontsize', 20,...
    'HorizontalAlignment', 'center')
text(ax(3),0.90,0.90,'R','Units','normalized', 'Fontsize', 20,...
    'HorizontalAlignment', 'center')

close(fig)
fig = gcf;

end
