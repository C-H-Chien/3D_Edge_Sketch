function finalPairIndx = getFinalHypothesisPair(Epipolar_Coeffs, edgels_HYPO2, max_support_indx, params)

    %> Code Description: 
    
    Edge_Pts = edgels_HYPO2(max_support_indx, 1:2);
    
    %> Vectorize the code
    Ap = Epipolar_Coeffs.A.*Edge_Pts(:,1);
    Bp = Epipolar_Coeffs.B.*Edge_Pts(:,2);
    numerOfDist = abs(Ap + Bp + Epipolar_Coeffs.C);
    denomOfDist = Epipolar_Coeffs.A.^2 + Epipolar_Coeffs.B.^2;
    denomOfDist = sqrt(denomOfDist);
    dist = numerOfDist./denomOfDist;
    
    [min_dist_val, min_dist_indx] = min(dist);
    if min_dist_val < params.SUPPORT_DIST_THRESH
        finalPairIndx = max_support_indx(min_dist_indx, 1);
    else
        finalPairIndx = 0;
    end
    
end