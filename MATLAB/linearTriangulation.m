function Gamma = linearTriangulation(N, pts_meters, Rs, Ts)
    %> Linear Triangulation Across Multiple (N) Views
    %
    %> Inputs:
    %    N:          Number of views
    %    pts_meters: Correspondences of points in meters of size 2xN; e.g.,
    %                for two views, pts_meters = [x_v1, x_v2;
    %                                             y_v1, y_v2]
    %    Rs:         Relative rotation matrix w.r.t. the first camera of 
    %                size 3x3x(N-1), i.e., Rs(:,:,1) is the relative 
    %                rotation matrix of view 2 w.r.t. view 1
    %    Ts:         Relative translation vector w.r.t. the first camera of 
    %                size 3x(N-1), i.e., Ts(:,:,1) is the relative 
    %                translation vector of view 2 w.r.t. view 1
    %
    %> Output:
    %    Gamma:      3D triangulated point located under the first camera
    %                coordinate
    %
    %> Author: Chiang-Heng Chien
    %> Last modified: April 10th, 2023
    
    %> Matrix A
    A = zeros(2*N, 4);
    
    %> We are computing 3D point under the first camera coordinate,
    %  so R1 is an identity matrix and T1 is all zeros
    A(1,1) = 0.0; A(1,2) = -1.0; A(1,3) = pts_meters(2,1); A(1,4) = 0.0;
    A(2,1) = 1.0; A(2,2) = 0.0; A(2,3)  = -pts_meters(1,1); A(2,4) = 0.0;
    
    row_cnter = 3;
    for p = 1:N-1
        %> Rotation and translation w.r.t. the first camera view
        Rp = Rs(:,:,p);
        Tp = Ts(:,p);
        mpx = pts_meters(1,p+1);
        mpy = pts_meters(2,p+1);
        
        r1 = Rp(1,1); r2 = Rp(1,2); r3 = Rp(1,3); t1 = Tp(1);
        r4 = Rp(2,1); r5 = Rp(2,2); r6 = Rp(2,3); t2 = Tp(2);
        r7 = Rp(3,1); r8 = Rp(3,2); r9 = Rp(3,3); t3 = Tp(3);
        
        A(row_cnter,   1) = mpy * r7 - r4;
        A(row_cnter,   2) = mpy * r8 - r5; 
        A(row_cnter,   3) = mpy * r9 - r6; 
        A(row_cnter,   4) = mpy * t3 - t2;
        A(row_cnter+1, 1) = r1 - mpx * r7; 
        A(row_cnter+1, 2) = r2 - mpx * r8; 
        A(row_cnter+1, 3) = r3 - mpx * r9; 
        A(row_cnter+1, 4) = t1 - mpx * t3;
        row_cnter = row_cnter + 2;
    end
    
    %> Triangulation
    ATA = A' * A;
    [~,~,V] = svd(ATA);
    p3d = V(:, end);
            
    %> Normalize X, Y, and Z by the last element
    p3d = p3d / p3d(end);
    Gamma = p3d(1:3,1);
end