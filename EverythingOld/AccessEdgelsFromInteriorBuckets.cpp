void getInteriorBuckets(
  const vgl_polygon_CH<double> &p, bool boundary_in, 
  std::vector<Eigen::Vector2d> &InteriorBucketCoordinates
)
{
  // iterate
  vgl_polygon_scan_iterator_CH<double> it(p, boundary_in); 

  //std::cout << "Interior points:\n";
  for (it.reset(); it.next(); ) {
    int y = it.scany();
    for (int x = it.startx(); x <= it.endx(); ++x) {
      //std::cout << "(" << x << "," << y << ") ";
      Eigen::Vector2d Bucket_Coordinate;
      Bucket_Coordinate << x, y;
      InteriorBucketCoordinates.push_back( Bucket_Coordinate );
    }
  }
  //std::cout << std::endl;
}

void getEdgelsFromInteriorQuadrilateral( 
  const subpixel_point_set &sp_pts, 
  const std::vector<Eigen::Vector2d> &InteriorBucketCoordinates,
  std::vector< unsigned > &Edgel_Indices 
)
{
  //> Traverse all buckets inside the quadrilateral
  //std::cout << "Number of interior bucket coordinates: " << InteriorBucketCoordinates.size()<<std::endl;
  //std::cout<<"bucket coordinates(starts from 0) are shown below: "<<std::endl;
  for (int bi = 0; bi < InteriorBucketCoordinates.size(); bi++) {
    
    unsigned const i_col = InteriorBucketCoordinates[bi](0);  //> x
    unsigned const i_row = InteriorBucketCoordinates[bi](1);  //> y

    //std::cout<< "coordinate " << bi << ": "<< i_col<< ", "<< i_row <<std::endl;
    //std::cout<< i_col << ", "<< i_row << ";" <<std::endl;

    //> Ignore if bucket coordinate exceeds image boundary
    if (i_row >= sp_pts.nrows() || i_col >= sp_pts.ncols()) continue;

    //> Traverse all edgels inside the bucket
    for (unsigned k = 0; k < sp_pts.cells()[i_row][i_col].size(); ++k) {
      unsigned const p2_idx = sp_pts.cells()[i_row][i_col][k];
      //std::cout<< "inlier edge index(starts from 0): " << p2_idx <<std::endl;
      Edgel_Indices.push_back(p2_idx);
    }
    //std::cout << "sp_pts.cells()[i_row][i_col].size(): " << sp_pts.cells()[i_row][i_col].size() <<std::endl;
  }
  //std::cout << "number of edges in this quadrilateral found by bucketing: "<< Edgel_Indices.size() <<std::endl;
}