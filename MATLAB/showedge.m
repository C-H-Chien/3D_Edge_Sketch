len = 3;
tgt1 = TO_Edges_Tangents{params.HYPO1_VIEW_INDX,1};
edgel_tgt1  = tgt1(edge_idx, 1:2);
hypo2_idx   = sort(hypo2_idx);
tgt2 = TO_Edges_Tangents{params.HYPO2_VIEW_INDX,1};
subplot(1,2,1);
hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
imshow(uint8(hypo_img1));
hold on;
plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
hold on;
x1 = edgel_HYPO1(1,1)+len*edgel_tgt1(1,1);
x2 = edgel_HYPO1(1,1)-len*edgel_tgt1(1,1);
y1 = edgel_HYPO1(1,2)+len*edgel_tgt1(1,2);
y2 = edgel_HYPO1(1,2)-len*edgel_tgt1(1,2);
hold on;
hypo31_all = edgel_HYPO1(1,1:2)' - epipole_pix_view1(1:2,:);
b31_all    = hypo31_all(2,:)./hypo31_all(1,:);
c31_all    = edgel_HYPO1(1,2) - b31_all * edgel_HYPO1(1,1);
y31_min1 = b31_all*1 + c31_all;
y31_max1 = b31_all*cols + c31_all;
hold on;
hold on;
plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'rx', 'MarkerSize',10, 'LineWidth', 1);
title 'hypothesis view 1'

subplot(1,2,2);
imshow(uint8(hypo_img2));
title 'hypothesis view 2'
hold on;
plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
hold on;
cols = size(hypo_img2, 2);
yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*cols) ./ Epipolar_Coeffs.B;

line([0, cols], [yMin, yMax], 'Color', 'y', 'LineWidth', 0.5);
hold on;
hypo31_all = edgels_HYPO2(:,1:2)' - epipole_pix_view2(1:2,:);
b31_all    = hypo31_all(2,:)./hypo31_all(1,:);
[wed_min,wed_minidx] = mink(b31_all,1);
[wed_max,wed_maxidx] = maxk(b31_all,1);
c31_all    = edgels_HYPO2(wed_minidx,2) - b31_all(1,wed_minidx) * edgels_HYPO2(wed_minidx,1);
y31_min1 = b31_all(1,wed_minidx)*1 + c31_all;
y31_max1 = b31_all(1,wed_minidx)*cols + c31_all;
hold on;
line([1, cols], [y31_min1, y31_max1], 'Color', 'y', 'LineWidth', 1);
c31_all    = edgels_HYPO2(wed_maxidx,2) - b31_all(1,wed_maxidx) * edgels_HYPO2(wed_maxidx,1);
y31_min1 = b31_all(1,wed_maxidx)*1 + c31_all;
y31_max1 = b31_all(1,wed_maxidx)*cols + c31_all;
hold on;
line([1, cols], [y31_min1, y31_max1], 'Color', 'y', 'LineWidth', 1);
plot(edgels_HYPO2(:, 1), edgels_HYPO2(:, 2), 'mx', 'MarkerSize',7, 'LineWidth', 1);
title 'hypothesis view 2'