load('ABC\ABC0006_GT.mat')
x_coorGT     = GT_all(:,1);
y_coorGT     = GT_all(:,2);
z_coorGT     = GT_all(:,3);
Curve_Pt_All = load('ABC\MCS_00000006_delta_theta_15.mat').recs;
Curve_Pt     = [];

for(idx = 1:size(Curve_Pt_All,2))
    Curve_Pt = [Curve_Pt;Curve_Pt_All{idx}];
end

Curve_Pt  = Curve_Pt';
dx1_curve = ones(1,size(x_coorGT,1)).*Curve_Pt(1,:)'-x_coorGT';
dy1_curve = ones(1,size(y_coorGT,1)).*Curve_Pt(2,:)'-y_coorGT';
dz1_curve = ones(1,size(z_coorGT,1)).*Curve_Pt(3,:)'-z_coorGT';

dist3Dall_curve    = sqrt(dx1_curve.^2 + dy1_curve.^2 + dz1_curve.^2);
[dist3D_curve,  ~] = mink(dist3Dall_curve,  1, 2);
[dist3D_curve1, ~] = mink(dist3Dall_curve', 1, 2);

cur_thr   = 0.02;
TP_curve  = find(dist3D_curve  < cur_thr);
TP_curve1 = find(dist3D_curve1 < cur_thr);

P_curve   = size(TP_curve, 1)  / size(dist3Dall_curve, 1);
R_curve   = size(TP_curve1, 1) / size(GT_all, 1);