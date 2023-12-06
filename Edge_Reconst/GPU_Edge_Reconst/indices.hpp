
#define Edgels_H1(i,j)          Edgels_H1[(i) * 3 + (j)]

#define All_Rot(i,j,k)          All_Rot[(i) * 3 + (j) + (k) * 9]
#define All_Calib(i,j,k)        All_Calib[(i) * 3 + (j) + (k) * 9]



//> Consistency Check
#define Converted_Matrix(i,j,k)     Converted_Matrix[(i) * 3 + (j) + (k) * 9]





//#define All_Transl(i,k)         All_Transl[(i) + (k) * 3];        //> (t) is useless in this case
