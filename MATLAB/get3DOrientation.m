function GammaTangent_3D = get3DOrientation(gamma_1,gamma_2, K, params, R21)
point1            = [gamma_1(1,1:2)'; 1];
point2            = [gamma_2(1,1:2)'; 1];
if(params.multiK == 1)
    K1 = K(:,:,params.HYPO1_VIEW_INDX);
    K2 = K(:,:,params.HYPO2_VIEW_INDX);
    invK1 = inv(K1);
    invK2 = inv(K2);

    gamma1            = invK1*point1;
    tangent1          = [cos(gamma_1(1,3)); ...
                         sin(gamma_1(1,3));1];
    pt1_tgt_to_pixels = point1(1:2,1) + tangent1(1:2,1);
    pt1_tgt_to_meters = invK1 * [pt1_tgt_to_pixels; 1];

    gamma2            = invK2*point2;
    tangent2          = [cos(gamma_2(1,3)); ...
                         sin(gamma_2(1,3));1];
    pt2_tgt_to_pixels = point2(1:2,1) + tangent2(1:2,1);
    pt2_tgt_to_meters = invK2 * [pt2_tgt_to_pixels; 1];

else
    invK = inv(K);

    gamma1            = invK*point1;
    tangent1          = [cos(gamma_1(1,3)); ...
                         sin(gamma_1(1,3));1];
    pt1_tgt_to_pixels = point1(1:2,1) + tangent1(1:2,1);
    pt1_tgt_to_meters = invK * [pt1_tgt_to_pixels; 1];

    gamma2            = invK*point2;
    tangent2          = [cos(gamma_2(1,3)); ...
                         sin(gamma_2(1,3));1];
    pt2_tgt_to_pixels = point2(1:2,1) + tangent2(1:2,1);
    pt2_tgt_to_meters = invK * [pt2_tgt_to_pixels; 1];
end

tgt1_meters       = pt1_tgt_to_meters - gamma1;
tgt2_meters       = pt2_tgt_to_meters - gamma2;

t1gamma1 = cross(tgt1_meters,gamma1);
t2gamma2 = cross(tgt2_meters,gamma2);
R21t2gamma2 = transpose(R21)*t2gamma2;
GammaTangent = cross(t1gamma1,R21t2gamma2);
GammaTangent_3D = GammaTangent./ norm(GammaTangent);

end