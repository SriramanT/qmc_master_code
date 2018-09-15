//
//  QDT.h
//
//
//  Created by Francisco Brito on 18/06/2018.
//

#ifndef QDT_h
#define QDT_h

template< int N >
class QDT
{
    Eigen::MatrixXd R;
    Eigen::MatrixXd T;
    Eigen::MatrixXd P;
public:
    Eigen::MatrixXd QR_and_getQ(Eigen::MatrixXd toDecompose);
    Eigen::MatrixXd getD();
    Eigen::MatrixXd getT();
    QDT() : R(N, N), T(N, N), P(N, N) {
    };
};

template< int N >
Eigen::MatrixXd QDT<N>::QR_and_getQ( Eigen::MatrixXd toDecompose )
{
    Eigen::ColPivHouseholderQR< Eigen::MatrixXd > colPivQrHH( toDecompose );
    R = colPivQrHH.matrixQR().template triangularView<Eigen::Upper>();
    P = colPivQrHH.colsPermutation();
    return colPivQrHH.householderQ();
}

template< int N >
Eigen::MatrixXd QDT<N>::getD()
{
    return ( R.diagonal() ).asDiagonal();
}

template< int N >
Eigen::MatrixXd QDT<N>::getT()
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

#endif /* QDT_h */
