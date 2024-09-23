load('ABC\ABC0006_GT.mat')
x_coorGT  = GT_all(:,1);
y_coorGT  = GT_all(:,2);
z_coorGT  = GT_all(:,3);
% Curve_all_theta005 = load('ABC\ABC0006_theta30theta1N4.mat').Gamma1_all;
% Curve_all_theta008 = load('ABC\ABC0006_theta30theta2N4.mat').Gamma1_all;
% Curve_all_theta01  = load('ABC\ABC0006_theta30theta5N4.mat').Gamma1_all;
% Curve_all_theta03  = load('ABC\ABC0006_theta30theta10N4.mat').Gamma1_all;
% Curve_all_theta06  = load('ABC\ABC0006_theta30theta15N4.mat').Gamma1_all;
% Curve_all_theta12  = load('ABC\ABC0006_theta30theta30N4.mat').Gamma1_all;

Curve_all_theta01_c  = load('ABC\MCS_00000077_delta_theta_01.mat').recs;
Curve_all_theta02_c  = load('ABC\MCS_00000077_delta_theta_02.mat').recs;
Curve_all_theta05_c  = load('ABC\MCS_00000077_delta_theta_05.mat').recs;
Curve_all_theta10_c  = load('ABC\MCS_00000077_delta_theta_10.mat').recs;
Curve_all_theta15_c  = load('ABC\MCS_00000077_delta_theta_15.mat').recs;
Curve_all_theta30_c  = load('ABC\MCS_00000077_delta_theta_30.mat').recs;


Curve_all_theta01  = [];%load('ABC\ABC0006_theta30theta1N4.mat').Gamma1_all;
Curve_all_theta02  = [];%load('ABC\ABC0006_theta30theta2N4.mat').Gamma1_all;
Curve_all_theta05  = [];%load('ABC\ABC0006_theta30theta5N4.mat').Gamma1_all;
Curve_all_theta10  = [];%load('ABC\ABC0006_theta30theta10N4.mat').Gamma1_all;
Curve_all_theta15  = [];%load('ABC\ABC0006_theta30theta15N4.mat').Gamma1_all;
Curve_all_theta30  = [];%load('ABC\ABC0006_theta30theta30N4.mat').Gamma1_all;

for(idx1 = 1:size(Curve_all_theta01_c,2))
    Curve_all_theta01 = [Curve_all_theta01;Curve_all_theta01_c{idx1}];
end

for(idx1 = 1:size(Curve_all_theta02_c,2))
    Curve_all_theta02 = [Curve_all_theta02;Curve_all_theta02_c{idx1}];
end

for(idx1 = 1:size(Curve_all_theta05_c,2))
    Curve_all_theta05 = [Curve_all_theta05;Curve_all_theta05_c{idx1}];
end

for(idx1 = 1:size(Curve_all_theta10_c,2))
    Curve_all_theta10 = [Curve_all_theta10;Curve_all_theta10_c{idx1}];
end

for(idx1 = 1:size(Curve_all_theta15_c,2))
    Curve_all_theta15 = [Curve_all_theta15;Curve_all_theta15_c{idx1}];
end

for(idx1 = 1:size(Curve_all_theta30_c,2))
    Curve_all_theta30 = [Curve_all_theta30;Curve_all_theta30_c{idx1}];
end

Curve_all_theta01 = Curve_all_theta01';%load('ABC\ABC0006_theta30theta1N4.mat').Gamma1_all;
Curve_all_theta02 = Curve_all_theta02';%load('ABC\ABC0006_theta30theta2N4.mat').Gamma1_all;
Curve_all_theta05  = Curve_all_theta05';%load('ABC\ABC0006_theta30theta5N4.mat').Gamma1_all;
Curve_all_theta10  = Curve_all_theta10';%load('ABC\ABC0006_theta30theta10N4.mat').Gamma1_all;
Curve_all_theta15  = Curve_all_theta15';%load('ABC\ABC0006_theta30theta15N4.mat').Gamma1_all;
Curve_all_theta30  = Curve_all_theta30';%load('ABC\ABC0006_theta30theta30N4.mat').Gamma1_all;

dx1_theta01 = ones(1,size(x_coorGT,1)).*Curve_all_theta01(1,:)'-x_coorGT';
dy1_theta01 = ones(1,size(y_coorGT,1)).*Curve_all_theta01(2,:)'-y_coorGT';
dz1_theta01 = ones(1,size(z_coorGT,1)).*Curve_all_theta01(3,:)'-z_coorGT';

dist3Dall1_theta01 = sqrt(dx1_theta01.^2 + dy1_theta01.^2 + dz1_theta01.^2);
[dist3D_Curve_theta01,~]=mink(dist3Dall1_theta01,1,2);

dx1_theta02 = ones(1,size(x_coorGT,1)).*Curve_all_theta02(1,:)'-x_coorGT';
dy1_theta02 = ones(1,size(y_coorGT,1)).*Curve_all_theta02(2,:)'-y_coorGT';
dz1_theta02 = ones(1,size(z_coorGT,1)).*Curve_all_theta02(3,:)'-z_coorGT';

dist3Dall1_theta02 = sqrt(dx1_theta02.^2 + dy1_theta02.^2 + dz1_theta02.^2);
[dist3D_Curve_theta02,~]=mink(dist3Dall1_theta02,1,2);

dx1_theta05 = ones(1,size(x_coorGT,1)).*Curve_all_theta05(1,:)'-x_coorGT';
dy1_theta05 = ones(1,size(y_coorGT,1)).*Curve_all_theta05(2,:)'-y_coorGT';
dz1_theta05 = ones(1,size(z_coorGT,1)).*Curve_all_theta05(3,:)'-z_coorGT';

dist3Dall1_theta05 = sqrt(dx1_theta05.^2 + dy1_theta05.^2 + dz1_theta05.^2);
[dist3D_Curve_theta05,~]=mink(dist3Dall1_theta05,1,2);

dx1_theta10 = ones(1,size(x_coorGT,1)).*Curve_all_theta10(1,:)'-x_coorGT';
dy1_theta10 = ones(1,size(y_coorGT,1)).*Curve_all_theta10(2,:)'-y_coorGT';
dz1_theta10 = ones(1,size(z_coorGT,1)).*Curve_all_theta10(3,:)'-z_coorGT';

dist3Dall1_theta10 = sqrt(dx1_theta10.^2 + dy1_theta10.^2 + dz1_theta10.^2);
[dist3D_Curve_theta10,~]=mink(dist3Dall1_theta10,1,2);

dx1_theta15 = ones(1,size(x_coorGT,1)).*Curve_all_theta15(1,:)'-x_coorGT';
dy1_theta15 = ones(1,size(y_coorGT,1)).*Curve_all_theta15(2,:)'-y_coorGT';
dz1_theta15 = ones(1,size(z_coorGT,1)).*Curve_all_theta15(3,:)'-z_coorGT';

dist3Dall1_theta15 = sqrt(dx1_theta15.^2 + dy1_theta15.^2 + dz1_theta15.^2);
[dist3D_Curve_theta15,~]=mink(dist3Dall1_theta15,1,2);

dx1_theta30 = ones(1,size(x_coorGT,1)).*Curve_all_theta30(1,:)'-x_coorGT';
dy1_theta30 = ones(1,size(y_coorGT,1)).*Curve_all_theta30(2,:)'-y_coorGT';
dz1_theta30 = ones(1,size(z_coorGT,1)).*Curve_all_theta30(3,:)'-z_coorGT';

dist3Dall1_theta30 = sqrt(dx1_theta30.^2 + dy1_theta30.^2 + dz1_theta30.^2);
[dist3D_Curve_theta30,~]=mink(dist3Dall1_theta30,1,2);

% dx = ones(1,size(x_coorGT,1)).*MCS_all(:,1)-x_coorGT';
% dy = ones(1,size(y_coorGT,1)).*MCS_all(:,2)-y_coorGT';
% dz = ones(1,size(z_coorGT,1)).*MCS_all(:,3)-z_coorGT';

% dist3Dall = sqrt(dx.^2 + dy.^2 + dz.^2);
% [dist3D_Curve,~]=mink(dist3Dall,1,2);

range_for_two  = [min([dist3D_Curve_theta30]) max([dist3D_Curve_theta30])];
sample_PR_dist = (range_for_two(1,2) - range_for_two(1,1))/10;
sample_PR      = [0.005,0.006,0.007,0.008,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1];

TP_curvet01  = [];
TP_curvet02  = [];
TP_curvet05  = [];
TP_curvet10  = [];
TP_curvet15  = [];
TP_curvet30  = [];


% for(sample_i = 1: size(sample_PR,2))
    cur_thr         = 0.02;%sample_PR(1,sample_i);
    TP_curvet01_cur = find(dist3D_Curve_theta01 <cur_thr);
    TP_curvet02_cur = find(dist3D_Curve_theta02 <cur_thr);
    TP_curvet05_cur  = find(dist3D_Curve_theta05  <cur_thr);
    TP_curvet10_cur  = find(dist3D_Curve_theta10  <cur_thr);
    TP_curvet15_cur  = find(dist3D_Curve_theta15  <cur_thr);
    TP_curvet30_cur  = find(dist3D_Curve_theta30  <cur_thr);
    TP_curvet01     = [TP_curvet01, size(TP_curvet01_cur,1)];
    TP_curvet02     = [TP_curvet02, size(TP_curvet02_cur,1)];
    TP_curvet05      = [TP_curvet05,  size(TP_curvet05_cur,1)];
    TP_curvet10      = [TP_curvet10,  size(TP_curvet10_cur,1)];
    TP_curvet15      = [TP_curvet15,  size(TP_curvet15_cur,1)];
    TP_curvet30      = [TP_curvet30,  size(TP_curvet30_cur,1)];
% end

P_curve = [];
R_curve = [];


P_curvet01 = TP_curvet01 /size(dist3D_Curve_theta01, 1);
P_curvet02 = TP_curvet02 /size(dist3D_Curve_theta02, 1);
P_curvet05  = TP_curvet05  /size(dist3D_Curve_theta05,  1);
P_curvet10  = TP_curvet10  /size(dist3D_Curve_theta10,  1);
P_curvet15  = TP_curvet15  /size(dist3D_Curve_theta15,  1);
P_curvet30  = TP_curvet30  /size(dist3D_Curve_theta30,  1);
P_curve     = [P_curvet01;P_curvet02;P_curvet05;P_curvet10;P_curvet15;P_curvet30];

R_curvet01 = TP_curvet01 /size(GT_all, 1);
R_curvet02 = TP_curvet02 /size(GT_all, 1);
R_curvet05  = TP_curvet05  /size(GT_all, 1);
R_curvet10  = TP_curvet10  /size(GT_all, 1);
R_curvet15  = TP_curvet15  /size(GT_all, 1);
R_curvet30  = TP_curvet30  /size(GT_all, 1);
R_curve     = [R_curvet01;R_curvet02;R_curvet05;R_curvet10;R_curvet15;R_curvet30];

% F_curve = 2*TP_curve/(size(GT_all,1)+size(dist3D_Curve,1));
% F_curve  = 2*TP_curve/(size(GT_all,1)+size(dist3D_Curve,1));

% figure
% plot(sample_PR, P_curve,'LineWidth',2)
% hold on
% plot(sample_PR, P_curve,'LineWidth',2)
% legend ({'3D Edge Sketch' '3D Curve Sketch'},'FontSize',15)
% xlabel ('threshold','FontSize',15)
% ylabel ('Precision','FontSize',15)
% figure
% plot(sample_PR, R_curve,'LineWidth',2)
% hold on
% plot(sample_PR, R_curve,'LineWidth',2)
% legend ({'3D Edge Sketch' '3D Curve Sketch'},'FontSize',15)
% xlabel ('threshold','FontSize',15)
% ylabel ('Recall','FontSize',15)
[P_curve';R_curve']

lineColor = [1,0,0; ...
             0,1,0; ...
             1,1,0; ...
             0,0,1; ...
             1,0,0; ...
             0,1,1; ...
             1,0,1; ...
             0,0,0; ...
             0,0.4470,0.7410; ...
             0.8500,0.3250,0.0980; ...
             0.9290,0.6940,0.1250; ...
             0.4940,0.1840,0.5560; ...
             0.4660,0.6740,0.1880; ...
             0.3010,0.7450,0.9330; ...
             0.6350,0.0780,0.1840];

% figure
% for(i = 1:1: size(sample_PR,2))
    hold on
p = plot(R_curve, P_curve,'LineWidth',2,'Color',lineColor(7,:));
% hold on
% plot(R_curve, P_curve,'LineWidth',2)
% legend ({'3D Edge Sketch' '3D Curve Sketch'},'FontSize',15)

p.DataTipTemplate.DataTipRows(1) = dataTipTextRow('Presicion',P_curve);
p.DataTipTemplate.DataTipRows(2) = dataTipTextRow('Recall',R_curve);
% p.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Δθ',[1,2,5,10,15,30]);
p.DataTipTemplate.DataTipRows(3) = dataTipTextRow('N support',[4,6,8,10,12,16]);
% p.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Δ',[0.05,0.08, 0.1,0.3,0.6,1.2]);
% end
xlabel ('Recall','FontSize',15)
ylabel ('Precision','FontSize',15)
% title  ('Precision-Recall Curve (changing Δθ)','FontSize',15)
% title  ('Precision-Recall Curve (changing N support)','FontSize',15)
title  ('Precision-Recall Curve (changing Δ)','FontSize',15)
% legend ({'D=0.005','D=0.006','D=0.007','D=0.008','D=0.01','D=0.02'},'FontSize',10)
% p.DataTipTemplate.DataTipRows(3) = dataTipTextRow('N support',[4,6,8,10,12,16]);
% p.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Δ',[0.05,0.08, 0.1,0.3,0.6,1.2]);
% text(R_curve(1,1),P_curve(1,1),'Δ=0.05','FontSize',15)
% text(R_curve(2,1),P_curve(2,1),'Δ=0.08','FontSize',15)
% text(R_curve(3,1),P_curve(3,1),'Δ=0.1','FontSize',15)
% text(R_curve(4,1),P_curve(4,1),'Δ=0.3','FontSize',15)
% text(R_curve(5,1),P_curve(5,1),'Δ=0.6','FontSize',15)
% text(R_curve(6,1),P_curve(6,1),'Δ=0.12','FontSize',15)
% legend ({'Δθ=15°, N support=4'},'FontSize',15)
% % text(R_curve(1,1),P_curve(1,1),'θ=1°','FontSize',15)
% % text(R_curve(2,1),P_curve(2,1),'θ=2°','FontSize',15)
% % text(R_curve(3,1),P_curve(3,1),'θ=5°','FontSize',15)
% % text(R_curve(4,1),P_curve(4,1),'θ=10°','FontSize',15)
% % text(R_curve(5,1),P_curve(5,1),'θ=15°','FontSize',15)
% % text(R_curve(6,1),P_curve(6,1),'θ=30°','FontSize',15)
% % legend ({'Δ=0.3, N support=4'},'FontSize',15)