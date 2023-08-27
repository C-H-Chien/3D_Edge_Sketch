#ifndef subpixel_point_set_h
#define subpixel_point_set_h
// 
//> Description: a class for efficient query of points within an image domain
//> Credit: Edited based on mw_subpixel_point_set from LEMSVXL by Ricardo Fabbri
// 
// Modifications
//    Chiang-Heng Chien  23-08-27    Built and made changes on top of mw_subpixel_point_set.h 
//                                   from LEMSVXL, under lemsvxl/contrib/rfabbri/mw/ in order 
//                                   to fit both the data structure and the need of Edge-Based 
//                                   Reconstruction.
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)

#include <vector>

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

class subpixel_point_set {
public:

    //> Constructor
    subpixel_point_set(const Eigen::MatrixXd &Edgels)
    : pt_(Edgels), nrows_(0), ncols_(0), NNN(1) {}
    
    //> Destructor
    ~subpixel_point_set() {};

    //> points represented by index into input point set vector.
    //  cells[row][col][edgel_idx]
    const std::vector< std::vector< std::vector<unsigned> > > & cells() const { return cells_; }

    //> Whether the image has bucketed or not
    bool is_bucketed() const { return cells_.size() != 0; }

    unsigned npts()  const { return pt_.rows(); }
    unsigned nrows() const { return nrows_; }
    unsigned ncols() const { return ncols_; }

    //> Build a bucketing grid
    void build_bucketing_grid( unsigned nr, unsigned nc ) 
    {
        nrows_ = nr; 
        ncols_ = nc; 

        cells_.clear();

        //> build a ncols x nrows cell grid
        cells_.resize(nrows_);
        for (unsigned ii=0; ii < nrows_; ii++) {
            cells_[ii].resize(ncols_);
        }

        for (unsigned ii=0; ii < nrows_; ii++) {
            for (unsigned jj=0; jj < ncols_; jj++)
                cells_[ii][jj].clear();
        }

        //> Push edgels into buckets
        for (unsigned id=0; id < pt_.rows(); ++id) 
        {
            //> Get the cell coordinates for this edgel
            unsigned gx = (unsigned)(pt_(id, 1) / NNN);
            unsigned gy = (unsigned)(pt_(id, 0) / NNN);

            //> if the edgel position is within the grid, assign it to the right grid point
            if ( gx < ncols_ && gy < nrows_) cells_[gy][gx].push_back(id);
        }

        //std::cout << "Done building bucketing grid!" << std::endl;
    }

protected:
  const Eigen::MatrixXd &pt_;
  std::vector< std::vector< std::vector<unsigned> > > cells_;
  unsigned       nrows_;     //> size of the bucketing grid
  unsigned       ncols_;     //> size of the bucketing grid
  unsigned const NNN;        //> may be used in future if coarser grid is desired
};


#endif