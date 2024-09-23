recons_coor1_all1 = load("ABC\ABC0006_testGamma1s.mat").Gamma1s;

recons_coor1_all2 = load("ABC\ABC0006_testGamma1s2.mat").Gamma1s;

recons_coor1_all3 = load("ABC\ABC0006_testGamma1s3.mat").Gamma1s;

recons_coor1_all4 = load("ABC\ABC0006_testGamma1s4.mat").Gamma1s;

recons_coor1_all5 = load("ABC\ABC0006_testGamma1s5.mat").Gamma1s;
figure;
plot3(recons_coor1_all1(1,:)*1, recons_coor1_all1(2,:)*1, recons_coor1_all1(3,:)*1, 'b.');
hold on;
plot3(recons_coor1_all2(1,:)*1, recons_coor1_all2(2,:)*1, recons_coor1_all2(3,:)*1, 'g.');
hold on;
plot3(recons_coor1_all3(1,:)*1, recons_coor1_all3(2,:)*1, recons_coor1_all3(3,:)*1, 'm.');
hold on;
plot3(recons_coor1_all4(1,:)*1, recons_coor1_all4(2,:)*1, recons_coor1_all4(3,:)*1, 'Color',[0.8500 0.3250 0.0980],'Marker','.','LineStyle','none');
hold on;
plot3(recons_coor1_all5(1,:)*1, recons_coor1_all5(2,:)*1, recons_coor1_all5(3,:)*1, 'Color',[0.3010 0.7450 0.9330],'Marker','.','LineStyle','none');
hold off;
axis equal;
view(3)