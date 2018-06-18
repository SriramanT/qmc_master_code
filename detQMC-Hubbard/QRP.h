//
//  QRP.h
//  
//
//  Created by Francisco Brito on 18/06/2018.
//

#ifndef QRP_h
#define QRP_h

template< int N >
class QRP
{
    Eigen::Matrix<double, N, N> R;
    Eigen::Matrix<double, N, N> T;
    Eigen::Matrix<double, N, N> P;
public:
    Eigen::Matrix<double, N, N> QR_and_getQ(Eigen::Matrix<double, N, N> toDecompose);
    Eigen::Matrix<double, N, N> getD();
    Eigen::Matrix<double, N, N> getT();
};

template< int N >
Eigen::Matrix<double, N, N> QRP<N>::QR_and_getQ( Eigen::Matrix<double, N, N> toDecompose )
{
    Eigen::ColPivHouseholderQR< Eigen::Matrix<double, N, N> > colPivQrHH( toDecompose );
    R = colPivQrHH.matrixQR().template triangularView<Eigen::Upper>();
    P = colPivQrHH.colsPermutation();
    return colPivQrHH.householderQ();
}

template< int N >
Eigen::Matrix<double, N, N> QRP<N>::getD()
{
    return ( R.diagonal() ).asDiagonal();
}

template< int N >
Eigen::Matrix<double, N, N> QRP<N>::getT()
{
    T = Eigen::Matrix<double, N, N>::Identity();
    for (int a = 0; a < N; a++)
    {
        for (int b = a + 1; b < N; b++) // unit upper triangular -> loop starts in i + 1
        {
            T(a, b) = R(a, b) / R(a, a);
        }
    }
    T = T * P.transpose();
    return T;
}

#endif /* QRP_h */
