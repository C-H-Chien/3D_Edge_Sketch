#include "linearTriangulationUtil.hpp"

Eigen::Vector3d linearTriangulation(const int N,
                                    const std::vector<Eigen::Vector2d> pts, 
                                    const std::vector<Eigen::Matrix3d> & Rs,
                                    const std::vector<Eigen::Vector3d> & Ts,
                                    const std::vector<double> & K) {
    
    Eigen::MatrixXd A(2*N, 4);
    Eigen::MatrixXd ATA(4, 4);
    Eigen::Vector4d GAMMA;
    Eigen::Vector3d pt3D; //> Returned variable

    //> Convert points in pixels to points in meters
    std::vector<Eigen::Vector2d> pts_meters;
    Eigen::Matrix3d K_;
    K_<<K[2], 0.,   K[0], 0.,  K[3], K[1], 0.,  0.,   1.;
    // std::cout << "K_: " << K_ <<std::endl;
    for (int p = 0; p < N; p++) {
        Eigen::Vector2d gamma;
        Eigen::Vector3d homo_p{pts[p](0), pts[p](1), 1.0};
        Eigen::Vector3d p_bar = K_.inverse() * homo_p;
        gamma(0) = p_bar(0);
        gamma(1) = p_bar(1);
        pts_meters.push_back(gamma);
    }

    // std::cout << "pts_meters[0]: " << pts_meters[0] <<std::endl;
    // std::cout << "pts_meters[1]: " << pts_meters[1] <<std::endl;

    //> We are computing GAMMA under the first camera coordinate,
    //  so R1 is an identity matrix and T1 is all zeros
    A(0,0) = 0.0; A(0,1) = -1.0; A(0,2) = pts_meters[0](1); A(0,3) = 0.0;
    A(1,0) = 1.0; A(1,1) = 0.0; A(1,2) = -pts_meters[0](0); A(1,3) = 0.0;

    int row_cnter = 2;
    for (int p = 0; p < N-1; p++) {
        Eigen::Matrix3d Rp = Rs[p];
        Eigen::Vector3d Tp = Ts[p];
        Eigen::Vector2d mp = pts_meters[p+1];

        // std::cout << "Rp: " << Rp <<std::endl;
        
        double r1 = Rp(0,0), r2 = Rp(0,1), r3 = Rp(0,2), t1 = Tp(0);
        double r4 = Rp(1,0), r5 = Rp(1,1), r6 = Rp(1,2), t2 = Tp(1);
        double r7 = Rp(2,0), r8 = Rp(2,1), r9 = Rp(2,2), t3 = Tp(2);
        // std::cout << "r1: " << r1 <<std::endl;
        // std::cout << "r2: " << r2 <<std::endl;
        // std::cout << "r3: " << r3 <<std::endl;
        // std::cout << "r4: " << r4 <<std::endl;
        // std::cout << "r5: " << r5 <<std::endl;
        // std::cout << "r6: " << r6 <<std::endl;
        // std::cout << "r7: " << r7 <<std::endl;
        // std::cout << "r8: " << r8 <<std::endl;
        // std::cout << "r9: " << r9 <<std::endl;
        // std::cout << "t1: " << t1 <<std::endl;
        // std::cout << "t2: " << t2 <<std::endl;
        // std::cout << "t3: " << t3 <<std::endl;

        A(row_cnter,   0) = mp(1) * r7 - r4;
        A(row_cnter,   1) = mp(1) * r8 - r5; 
        A(row_cnter,   2) = mp(1) * r9 - r6; 
        A(row_cnter,   3) = mp(1) * t3 - t2;
        A(row_cnter+1, 0) = r1 - mp(0) * r7; 
        A(row_cnter+1, 1) = r2 - mp(0) * r8; 
        A(row_cnter+1, 2) = r3 - mp(0) * r9; 
        A(row_cnter+1, 3) = t1 - mp(0) * t3;
        row_cnter += 2;
    }

    // std::cout << "A: " << A <<std::endl;

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