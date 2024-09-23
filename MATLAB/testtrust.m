rangeplot = [203]%   204   208   214   215   216   221   222   229   230   238   239   247   255];
figure;
plot3(recons_coor2(:,1), -recons_coor2(:,2), -recons_coor2(:,3), 'b.', 'MarkerSize', 3, 'LineWidth', 1);
hold on;
plot3(recons_coor2(rangeplot,1), -recons_coor2(rangeplot,2), -recons_coor2(rangeplot,3), 'gx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
% title (['3D reconstruction result (accumulative, wedge = ± ', num2str(params.ANGLE_FOR_EPIPOLE),'°)']) 
title (['3D reconstruction result (single round after BA)'])