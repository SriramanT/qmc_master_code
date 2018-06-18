//
//  SVD.h
//  
//
//  Created by Francisco Brito on 18/06/2018.
//

#ifndef SVD_h
#define SVD_h

template< int N >
class SVD
{
    Eigen::Matrix<double, N, N> U;
    Eigen::Matrix<double, N, N> S;
    Eigen::Matrix<double, N, N> V;
public:
    void doSVD( Eigen::Matrix<double, N, N> toDecompose );
    Eigen::Matrix<double, N, N> getU();
    Eigen::Matrix<double, N, N> getS();
    Eigen::Matrix<double, N, N> getV();
};

template< int N >
void SVD<N>::doSVD( Eigen::Matrix<double, N, N> toDecompose )
{
    Eigen::JacobiSVD< Eigen::Matrix<double, N, N> > svd( toDecompose , Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    S = ( svd.singularValues() ).asDiagonal();
    V = svd.matrixV().transpose();
}

template< int N >
Eigen::Matrix<double, N, N> SVD<N>::getU()
{
    return U;
}

template< int N >
Eigen::Matrix<double, N, N> SVD<N>::getS()
{
    return S;
}

template< int N >
Eigen::Matrix<double, N, N> SVD<N>::getV()
{
    return V;
}

#endif /* SVD_h */
