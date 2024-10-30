//> Dependencies: 
//  - Eigen 3.3.2
//  - OpenCV 3.4.1 or higher
//
//> Descriptions: linear triangulation of point correspondence across N views
//> Inputs:
//   - N:          Number of views
//   - pts_meters: Correspondences of points in pixels in the form of cv::Point2d and stacked by std::vector structure.
//                 For example: cv::Point2d pt2D( x, y );
//                              std::vector pts.push_back(pt2D);
//   - Rs:         Relative rotation matrix w.r.t. the first camera in the form of Eigen::Matrix3d and stacked by std::vector structure.
//                 The first one, i.e., Rs[0], is the relative rotation matrix of the second view w.r.t. the first view.
//   - Ts:         Relative translation vector w.r.t. the first camera in the form of Eigen::Vector3d and stacked by std::vector structure 
//                 The first one, i.e., Ts[0], is the relative translation vector of the second view w.r.t. the first view.
//   - K:          Camera intrinsic matrix. Elements inserted in a std::vector<double> structure.
//                 Specifically, K[0] = cx, K[1] = cy, K[2] = fx, K[3] = fy, assuming the skew is zero.
//                 cx and cy are the optic center, while fx and fy are the focal lengths.
//
//> Output: 3D triangulated point 'pt3D' located under the first camera coordinate, in the form of Eigen::Vector3d.
//
//> Author: Chiang-Heng Chien
//> Last modified: Feb. 20th, 2023

Eigen::Vector3d linearTriangulation(const int N,
                                    const std::vector<cv::Point2d> pts, 
                                    const std::vector<Eigen::Matrix3d> & Rs,
                                    const std::vector<Eigen::Vector3d> & Ts,
                                    const std::vector<double> & K) {
    
    Eigen::MatrixXd A(2*N, 4);
    Eigen::MatrixXd ATA(4, 4);
    Eigen::Vector4d GAMMA;
    Eigen::Vector3d pt3D; //> Returned variable

    //> Convert points in pixels to points in meters
    std::vector<cv::Point2d> pts_meters;
    Eigen::Matrix3d K_{{K[2], 0.,   K[0]},
                       {0.,  K[3], K[1]},
                       {0.,  0.,   1.}};
    for (int p = 0; p < N; p++) {
        cv::Point2d gamma;
        Eigen::Vector3d homo_p{pts[p].x, pts[p].y, 1.0};
        Eigen::Vector3d p_bar = K_.inverse() * homo_p;
        gamma.x = p_bar(0);
        gamma.y = p_bar(1);
        pts_meters.push_back(gamma);
    }

    //> We are computing GAMMA under the first camera coordinate,
    //  so R1 is an identity matrix and T1 is all zeros
    A(0,0) = 0.0; A(0,1) = -1.0; A(0,2) = pts_meters[0].y; A(0,3) = 0.0;
    A(1,0) = 1.0; A(1,1) = 0.0; A(1,2) = -pts_meters[0].x; A(1,3) = 0.0;

    int row_cnter = 2;
    for (int p = 0; p < N-1; p++) {
        Eigen::Matrix3d Rp = Rs[p];
        Eigen::Vector3d Tp = Ts[p];
        cv::Poits_meters[p+1];
        
        double r1 = Rp(0,0), r2 = Rp(0,1), r3 = Rp(0,2), t1 = Tp(0);
        double r4 = Rp(1,0), r5 = Rp(1,1), r6 = Rp(1,2), t2 = Tp(1);
        double r7 = Rp(2,0), r8 = Rp(2,1), r9 = Rp(2,2), t3 = Tp(2);

        A(row_cnter,   0) = mp.y * r7 - r4;
        A(row_cnter,   1) = mp.y * r8 - r5; 
        A(row_cnter,   2) = mp.y * r9 - r6; 
        A(row_cnter,   3) = mp.y * t3 - t2;
        A(row_cnter+1, 0) = r1 - mp.x * r7; 
        A(row_cnter+1, 1) = r2 - mp.x * r8; 
        A(row_cnter+1, 2) = r3 - mp.x * r9; 
        A(row_cnter+1, 3) = t1 - mp.x * t3;
        row_cnter += 2;
    }

    //> Solving the homogeneous linear system and divide the first three rows with the last element
    ATA = A.transpose() * A;
    GAMMA = ATA.jacobiSvd(Eigen::ComputeFullV).matrixV().col( ATA.rows() - 1 );
    GAMMA[0] /= GAMMA[3];
    GAMMA[1] /= GAMMA[3];
    GAMMA[2] /= GAMMA[3];
    
    //> Assign GAMMA to the returned point
    pt3D[0] = GAMMA[0];
    pt3D[1] = GAMMA[1];
    pt3D[2] = GAMMA[2];
    
    return pt3D;

    /*
    template<typename matrix_t, typename vector_t>
    void solveNullspaceLU(const matrix_t& A, vector_t& x){
        x = A.fullPivLu().kernel();
        x.normalize();
    }

    template<typename matrix_t, typename vector_t>
    void solveNullspaceQR(const matrix_t& A, vector_t& x){
        auto qr = A.transpose().colPivHouseholderQr();
        matrix_t Q = qr.householderQ();
        x = Q.col(A.rows() - 1);
        x.normalize();
    }

    template<typename matrix_t, typename vector_t>
    void solveNullspaceSVD(const matrix_t& A, vector_t& x){
        x = A.jacobiSvd(Eigen::ComputeFullV).matrixV().col( A.rows() - 1 );
        x.normalize();
    }
    */
}
