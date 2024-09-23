function [edgels_HYPO2_corrected] = edgelsHYPO2correct_1(Epipolar_Coeffs, edgels_HYPO2)

edgels_HYPO2_corrected = edgels_HYPO2;

% match1 = [img1_points, ones(size(img1_points, 1), 1)]';
% match2 = [img2_points, ones(size(img1_points, 1), 1)]';

% InitialReprojectionErrors = zeros(size(edgels_HYPO2_corrected, 1), 1);
% OptimizedReprojectionErrors = zeros(size(edgels_HYPO2_corrected, 1), 1);

a_line = -Epipolar_Coeffs.A./Epipolar_Coeffs.B;
b_line = -1;
c_line = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;

for i = 1:size(edgels_HYPO2_corrected, 1)
    hypo2_current = edgels_HYPO2(i,:);
    slop_current  = tan(hypo2_current(1,3));
    tgt_current   = [cos(hypo2_current(1,3)), ...
                     sin(hypo2_current(1,3))];
%     x1 = hypo2_current(1,1)+len*tgt_current(1,1);
%     x2 = hypo2_current(1,1)-len*tgt_current(1,1);
%     y1 = hypo2_current(1,2)+len*tgt_current(1,2);
%     y2 = hypo2_current(1,2)-len*tgt_current(1,2);

    a_current = slop_current;
    b_current = -1;
    c_current = -(slop_current*hypo2_current(1,1)-hypo2_current(1,2));

    correctedpt = [(b_line*c_current-b_current*c_line)/(a_line*b_current-a_current*b_line), (c_line*a_current-c_current*a_line)/(a_line*b_current-a_current*b_line)];
    edgels_HYPO2_corrected(i, 1:2) = correctedpt;
end

end