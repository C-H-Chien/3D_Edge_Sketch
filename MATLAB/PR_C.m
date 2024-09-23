load('ABC\ABC0006_GT.mat')
x_coorGT  = GT_all(:,1);
y_coorGT  = GT_all(:,2);
z_coorGT  = GT_all(:,3);
Edge_all_delta005 = load('ABC\ABC0006_testGamma1sall.mat').Gamma1_all;
Edge_all_delta008 = load('ABC\ABC0006_delta03theta2N2.mat').Gamma1_all;
Edge_all_delta01  = load('ABC\ABC0006_delta03theta5N2.mat').Gamma1_all;
Edge_all_delta03  = load('ABC\ABC0006_delta03theta10N2.mat').Gamma1_all;
Edge_all_delta06  = load('ABC\ABC0006_delta03theta15N4_4vp.mat').Gamma1_all;
Edge_all_delta12  = load('ABC\ABC0006_delta03theta30N2.mat').Gamma1_all;

% % % % % Edge_all_delta005_c = load('ABC\MCS_00000006_delta_theta_01.mat').recs;
% % % % % Edge_all_delta008_c = load('ABC\MCS_00000006_delta_theta_02.mat').recs;
% % % % % Edge_all_delta01_c  = load('ABC\MCS_00000006_delta_theta_05.mat').recs;
% % % % % Edge_all_delta03_c  = load('ABC\MCS_00000006_delta_theta_10.mat').recs;
% % % % % Edge_all_delta06_c  = load('ABC\MCS_00000006_delta_theta_15.mat').recs;
% % % % % Edge_all_delta12_c  = load('ABC\MCS_00000006_delta_theta_30.mat').recs;
% % % % % 
% % % % % 
% % % % % Edge_all_delta005 = [];%load('ABC\ABC0006_delta12theta1N4.mat').Gamma1_all;
% % % % % Edge_all_delta008 = [];%load('ABC\ABC0006_delta12theta2N4.mat').Gamma1_all;
% % % % % Edge_all_delta01  = [];%load('ABC\ABC0006_delta12theta5N4.mat').Gamma1_all;
% % % % % Edge_all_delta03  = [];%load('ABC\ABC0006_delta12theta10N4.mat').Gamma1_all;
% % % % % Edge_all_delta06  = [];%load('ABC\ABC0006_delta12theta15N4.mat').Gamma1_all;
% % % % % Edge_all_delta12  = [];%load('ABC\ABC0006_delta12theta30N4.mat').Gamma1_all;
% % % % % 
% % % % % for(idx1 = 1:size(Edge_all_delta005_c,2))
% % % % %     Edge_all_delta005 = [Edge_all_delta005;Edge_all_delta005_c{idx1}];
% % % % % end
% % % % % 
% % % % % for(idx1 = 1:size(Edge_all_delta008_c,2))
% % % % %     Edge_all_delta008 = [Edge_all_delta008;Edge_all_delta008_c{idx1}];
% % % % % end
% % % % % 
% % % % % for(idx1 = 1:size(Edge_all_delta01_c,2))
% % % % %     Edge_all_delta01 = [Edge_all_delta01;Edge_all_delta01_c{idx1}];
% % % % % end
% % % % % 
% % % % % for(idx1 = 1:size(Edge_all_delta03_c,2))
% % % % %     Edge_all_delta03 = [Edge_all_delta03;Edge_all_delta03_c{idx1}];
% % % % % end
% % % % % 
% % % % % for(idx1 = 1:size(Edge_all_delta06_c,2))
% % % % %     Edge_all_delta06 = [Edge_all_delta06;Edge_all_delta06_c{idx1}];
% % % % % end
% % % % % 
% % % % % for(idx1 = 1:size(Edge_all_delta12_c,2))
% % % % %     Edge_all_delta12 = [Edge_all_delta12;Edge_all_delta12_c{idx1}];
% % % % % end
% % % % % 
% % % % % Edge_all_delta005 = Edge_all_delta005';%load('ABC\ABC0006_delta12theta1N4.mat').Gamma1_all;
% % % % % Edge_all_delta008 = Edge_all_delta008';%load('ABC\ABC0006_delta12theta2N4.mat').Gamma1_all;
% % % % % Edge_all_delta01  = Edge_all_delta01';%load('ABC\ABC0006_delta12theta5N4.mat').Gamma1_all;
% % % % % Edge_all_delta03  = Edge_all_delta03';%load('ABC\ABC0006_delta12theta10N4.mat').Gamma1_all;
% % % % % Edge_all_delta06  = Edge_all_delta06';%load('ABC\ABC0006_delta12theta15N4.mat').Gamma1_all;
% % % % % Edge_all_delta12  = Edge_all_delta12';%load('ABC\ABC0006_delta12theta30N4.mat').Gamma1_all;

dx1_delta005 = ones(1,size(x_coorGT,1)).*Edge_all_delta005(1,:)'-x_coorGT';
dy1_delta005 = ones(1,size(y_coorGT,1)).*Edge_all_delta005(2,:)'-y_coorGT';
dz1_delta005 = ones(1,size(z_coorGT,1)).*Edge_all_delta005(3,:)'-z_coorGT';

dist3Dall1_delta005 = sqrt(dx1_delta005.^2 + dy1_delta005.^2 + dz1_delta005.^2);
[dist3D_Edge_delta005,~]=mink(dist3Dall1_delta005,1,2);

dx1_delta008 = ones(1,size(x_coorGT,1)).*Edge_all_delta008(1,:)'-x_coorGT';
dy1_delta008 = ones(1,size(y_coorGT,1)).*Edge_all_delta008(2,:)'-y_coorGT';
dz1_delta008 = ones(1,size(z_coorGT,1)).*Edge_all_delta008(3,:)'-z_coorGT';

dist3Dall1_delta008 = sqrt(dx1_delta008.^2 + dy1_delta008.^2 + dz1_delta008.^2);
[dist3D_Edge_delta008,~]=mink(dist3Dall1_delta008,1,2);

dx1_delta01 = ones(1,size(x_coorGT,1)).*Edge_all_delta01(1,:)'-x_coorGT';
dy1_delta01 = ones(1,size(y_coorGT,1)).*Edge_all_delta01(2,:)'-y_coorGT';
dz1_delta01 = ones(1,size(z_coorGT,1)).*Edge_all_delta01(3,:)'-z_coorGT';

dist3Dall1_delta01 = sqrt(dx1_delta01.^2 + dy1_delta01.^2 + dz1_delta01.^2);
[dist3D_Edge_delta01,~]=mink(dist3Dall1_delta01,1,2);

dx1_delta03 = ones(1,size(x_coorGT,1)).*Edge_all_delta03(1,:)'-x_coorGT';
dy1_delta03 = ones(1,size(y_coorGT,1)).*Edge_all_delta03(2,:)'-y_coorGT';
dz1_delta03 = ones(1,size(z_coorGT,1)).*Edge_all_delta03(3,:)'-z_coorGT';

dist3Dall1_delta03 = sqrt(dx1_delta03.^2 + dy1_delta03.^2 + dz1_delta03.^2);
[dist3D_Edge_delta03,~]=mink(dist3Dall1_delta03,1,2);

dx1_delta06 = ones(1,size(x_coorGT,1)).*Edge_all_delta06(1,:)'-x_coorGT';
dy1_delta06 = ones(1,size(y_coorGT,1)).*Edge_all_delta06(2,:)'-y_coorGT';
dz1_delta06 = ones(1,size(z_coorGT,1)).*Edge_all_delta06(3,:)'-z_coorGT';

dist3Dall1_delta06 = sqrt(dx1_delta06.^2 + dy1_delta06.^2 + dz1_delta06.^2);
[dist3D_Edge_delta06,~]=mink(dist3Dall1_delta06,1,2);

dx1_delta12 = ones(1,size(x_coorGT,1)).*Edge_all_delta12(1,:)'-x_coorGT';
dy1_delta12 = ones(1,size(y_coorGT,1)).*Edge_all_delta12(2,:)'-y_coorGT';
dz1_delta12 = ones(1,size(z_coorGT,1)).*Edge_all_delta12(3,:)'-z_coorGT';

dist3Dall1_delta12 = sqrt(dx1_delta12.^2 + dy1_delta12.^2 + dz1_delta12.^2);
[dist3D_Edge_delta12,~]=mink(dist3Dall1_delta12,1,2);

% dx = ones(1,size(x_coorGT,1)).*MCS_all(:,1)-x_coorGT';
% dy = ones(1,size(y_coorGT,1)).*MCS_all(:,2)-y_coorGT';
% dz = ones(1,size(z_coorGT,1)).*MCS_all(:,3)-z_coorGT';

% dist3Dall = sqrt(dx.^2 + dy.^2 + dz.^2);
% [dist3D_Curve,~]=mink(dist3Dall,1,2);

range_for_two  = [min([dist3D_Edge_delta12]) max([dist3D_Edge_delta12])];
sample_PR_dist = (range_for_two(1,2) - range_for_two(1,1))/10;
sample_PR      = [0.005,0.006,0.007,0.008,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1];

TP_edged005  = [];
TP_edged008  = [];
TP_edged01  = [];
TP_edged03  = [];
TP_edged06  = [];
TP_edged12  = [];


% for(sample_i = 1: size(sample_PR,2))
    cur_thr         = 0.02;%sample_PR(1,sample_i);
    TP_edged005_cur = find(dist3D_Edge_delta005 <cur_thr);
    TP_edged008_cur = find(dist3D_Edge_delta008 <cur_thr);
    TP_edged01_cur  = find(dist3D_Edge_delta01  <cur_thr);
    TP_edged03_cur  = find(dist3D_Edge_delta03  <cur_thr);
    TP_edged06_cur  = find(dist3D_Edge_delta06  <cur_thr);
    TP_edged12_cur  = find(dist3D_Edge_delta12  <cur_thr);
    TP_edged005     = [TP_edged005, size(TP_edged005_cur,1)];
    TP_edged008     = [TP_edged008, size(TP_edged008_cur,1)];
    TP_edged01      = [TP_edged01,  size(TP_edged01_cur,1)];
    TP_edged03      = [TP_edged03,  size(TP_edged03_cur,1)];
    TP_edged06      = [TP_edged06,  size(TP_edged06_cur,1)];
    TP_edged12      = [TP_edged12,  size(TP_edged12_cur,1)];
% end

P_edge = [];
R_edge = [];


P_edged005 = TP_edged005 /size(dist3D_Edge_delta005, 1);
P_edged008 = TP_edged008 /size(dist3D_Edge_delta008, 1);
P_edged01  = TP_edged01  /size(dist3D_Edge_delta01,  1);
P_edged03  = TP_edged03  /size(dist3D_Edge_delta03,  1);
P_edged06  = TP_edged06  /size(dist3D_Edge_delta06,  1);
P_edged12  = TP_edged12  /size(dist3D_Edge_delta12,  1);
P_edge     = [P_edged005;P_edged008;P_edged01;P_edged03;P_edged06;P_edged12];

R_edged005 = TP_edged005 /size(GT_all, 1);
R_edged008 = TP_edged008 /size(GT_all, 1);
R_edged01  = TP_edged01  /size(GT_all, 1);
R_edged03  = TP_edged03  /size(GT_all, 1);
R_edged06  = TP_edged06  /size(GT_all, 1);
R_edged12  = TP_edged12  /size(GT_all, 1);
R_edge     = [R_edged005;R_edged008;R_edged01;R_edged03;R_edged06;R_edged12];

% F_curve = 2*TP_curve/(size(GT_all,1)+size(dist3D_Curve,1));
% F_edge  = 2*TP_edge/(size(GT_all,1)+size(dist3D_Curve,1));

% figure
% plot(sample_PR, P_edge,'LineWidth',2)
% hold on
% plot(sample_PR, P_curve,'LineWidth',2)
% legend ({'3D Edge Sketch' '3D Curve Sketch'},'FontSize',15)
% xlabel ('threshold','FontSize',15)
% ylabel ('Precision','FontSize',15)
% figure
% plot(sample_PR, R_edge,'LineWidth',2)
% hold on
% plot(sample_PR, R_curve,'LineWidth',2)
% legend ({'3D Edge Sketch' '3D Curve Sketch'},'FontSize',15)
% xlabel ('threshold','FontSize',15)
% ylabel ('Recall','FontSize',15)
[P_edge';R_edge']

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
p = plot(R_edge, P_edge,'LineWidth',2,'Color',lineColor(7,:));
% hold on
% plot(R_curve, P_curve,'LineWidth',2)
% legend ({'3D Edge Sketch' '3D Curve Sketch'},'FontSize',15)

p.DataTipTemplate.DataTipRows(1) = dataTipTextRow('Presicion',P_edge);
p.DataTipTemplate.DataTipRows(2) = dataTipTextRow('Recall',R_edge);
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
% text(R_edge(1,1),P_edge(1,1),'Δ=0.05','FontSize',15)
% text(R_edge(2,1),P_edge(2,1),'Δ=0.08','FontSize',15)
% text(R_edge(3,1),P_edge(3,1),'Δ=0.1','FontSize',15)
% text(R_edge(4,1),P_edge(4,1),'Δ=0.3','FontSize',15)
% text(R_edge(5,1),P_edge(5,1),'Δ=0.6','FontSize',15)
% text(R_edge(6,1),P_edge(6,1),'Δ=0.12','FontSize',15)
% legend ({'Δθ=15°, N support=4'},'FontSize',15)
% % text(R_edge(1,1),P_edge(1,1),'θ=1°','FontSize',15)
% % text(R_edge(2,1),P_edge(2,1),'θ=2°','FontSize',15)
% % text(R_edge(3,1),P_edge(3,1),'θ=5°','FontSize',15)
% % text(R_edge(4,1),P_edge(4,1),'θ=10°','FontSize',15)
% % text(R_edge(5,1),P_edge(5,1),'θ=15°','FontSize',15)
% % text(R_edge(6,1),P_edge(6,1),'θ=30°','FontSize',15)
% % legend ({'Δ=0.3, N support=4'},'FontSize',15)