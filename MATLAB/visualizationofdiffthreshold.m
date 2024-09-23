load('ABC0006_3D.mat')
figure(1)
subplot(2,2,1)
p3 = scatter3(Edge_all(1,:),       Edge_all(2,:),       Edge_all(3,:), ...
    1,'MarkerEdgeColor','r','Marker','x');
axis equal
grid off
xlim([min(Edge_all_nore(1,:)) max(Edge_all_nore(1,:))])
ylim([min(Edge_all_nore(2,:)) max(Edge_all_nore(2,:))])
zlim([min(Edge_all_nore(3,:)) max(Edge_all_nore(3,:))])
title ('3D Edge Sketch with delta=0.3, orentation threshold= 15°','FontSize',12)
subplot(2,2,2)
p3 = scatter3(Edge_all_nore(1,:),       Edge_all_nore(2,:),       Edge_all_nore(3,:), ...
    1,'MarkerEdgeColor','#0072BD','Marker','x');
axis equal
grid off
xlim([min(Edge_all_nore(1,:)) max(Edge_all_nore(1,:))])
ylim([min(Edge_all_nore(2,:)) max(Edge_all_nore(2,:))])
zlim([min(Edge_all_nore(3,:)) max(Edge_all_nore(3,:))])
title ('3D Edge Sketch with delta=0.3, orentation threshold= 180°','FontSize',12)
subplot(2,2,3)
p3 = scatter3(Edge_all_delta01(1,:),       Edge_all_delta01(2,:),       Edge_all_delta01(3,:), ...
    1,'MarkerEdgeColor','#D95319','Marker','x');
axis equal
grid off
xlim([min(Edge_all_nore(1,:)) max(Edge_all_nore(1,:))])
ylim([min(Edge_all_nore(2,:)) max(Edge_all_nore(2,:))])
zlim([min(Edge_all_nore(3,:)) max(Edge_all_nore(3,:))])
title ('3D Edge Sketch with with delta=0.1, orentation threshold= 15°','FontSize',12)
subplot(2,2,4)
p3 = scatter3(Edge_all_delta1(1,:),       Edge_all_delta1(2,:),       Edge_all_delta1(3,:), ...
    1,'MarkerEdgeColor','#7E2F8E','Marker','x');
axis equal
grid off
xlim([min(Edge_all_nore(1,:)) max(Edge_all_nore(1,:))])
ylim([min(Edge_all_nore(2,:)) max(Edge_all_nore(2,:))])
zlim([min(Edge_all_nore(3,:)) max(Edge_all_nore(3,:))])
title ('3D Edge Sketch with delta=1, orentation threshold= 15°','FontSize',12)