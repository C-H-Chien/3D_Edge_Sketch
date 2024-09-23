sample_PR = 0.01;
TP_edge  = [];
TP_curve = [];
TP_edged01  = [];
TP_edged05  = [];
TP_edged1  = [];

for(sample_i = 1: size(sample_PR,2))
    cur_thr      = sample_PR(1,sample_i);
    TP_edge_cur  = find(dist3D_Edge <cur_thr);
    TP_curve_cur = find(dist3D_Curve <cur_thr);
    TP_edged01_cur  = find(dist3D_Edge_delta01 <cur_thr);
    TP_edged05_cur = find(dist3D_Edge_delta05 <cur_thr);
    TP_edged1_cur  = find(dist3D_Edge_delta1 <cur_thr);
    TP_edge      = [TP_edge,    size(TP_edge_cur,1)];
    TP_curve     = [TP_curve,   size(TP_curve_cur,1)];
    TP_edged01   = [TP_edged01, size(TP_edged01_cur,1)];
    TP_edged05   = [TP_edged05, size(TP_edged05_cur,1)];
    TP_edged1    = [TP_edged1,  size(TP_edged1_cur,1)];
end

P_curve    = TP_curve/size(dist3D_Curve,1);
P_edge     = TP_edge/size(dist3D_Edge,1);
P_edged01  = TP_edged01/size(dist3D_Edge_delta01,1);
P_edged05  = TP_edged05/size(dist3D_Edge_delta05,1);
P_edged1   = TP_edged1/size(dist3D_Edge_delta1,1);

R_curve = TP_curve/size(GT_all,1);
R_edge  = TP_edge/size(GT_all,1);
R_edged01  = TP_edged01/size(GT_all,1);
R_edged05  = TP_edged05/size(GT_all,1);
R_edged1   = TP_edged1/size(GT_all,1);

R_all = [R_edged01;R_edge;R_edged05;R_edged1];
P_all = [P_edged01;P_edge;P_edged05;P_edged1];

figure;
plot(R_all, P_all,'LineWidth',2,"Color",[0 0.447 0.741])
xlabel ('Recall','FontSize',15)
ylabel ('Precision','FontSize',15)
title  ('Precision-Recall Curve','FontSize',15)