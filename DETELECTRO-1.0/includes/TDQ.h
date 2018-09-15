//
//  TDQ.h
//
//
//  Created by Francisco Brito on 18/06/2018.
//

#ifndef TDQ_h
#define TDQ_h

template< int N >
class TDQ
{
    Eigen::MatrixXd Q;
    Eigen::MatrixXd R;
    Eigen::MatrixXd T;
    Eigen::MatrixXd P;
public:
    Eigen::MatrixXd QR_and_getQ(Eigen::MatrixXd toDecompose);
    Eigen::MatrixXd getD();
    Eigen::MatrixXd getT();
    TDQ() : Q(N, N), R(N, N), T(N, N), P(N, N)  {
    };
};

template< int N >
Eigen::MatrixXd TDQ<N>::QR_and_getQ( Eigen::MatrixXd toDecompose )
{
    Eigen::ColPivHouseholderQR< Eigen::MatrixXd > colPivQrHH( toDecompose.colwise().reverse().transpose() );
    Q = colPivQrHH.householderQ();
    R = colPivQrHH.matrixQR().template triangularView<Eigen::Upper>();
    R = R.transpose().colwise().reverse().rowwise().reverse().eval();
    P = colPivQrHH.colsPermutation();
    P = P.colwise().reverse().rowwise().reverse().eval();
    return Q.transpose().colwise().reverse();
}

template< int N >
Eigen::MatrixXd TDQ<N>::getD()
{
    return ( R.diagonal() ).asDiagonal();
}

template< int N >
Eigen::MatrixXd TDQ<N>::getT()
{
    T = Eigen::Matrix<double, N, N>::Identity();
    for (int a = 0; a < N; a++)
    {
        for (int b = a + 1; b < N; b++) // unit upper triangular -> loop starts in i + 1
        {
            T(a, b) = R(a, b) / R(b, b);
        }
    }
    T = P * T;
    return T;
}

#endif /* TDQ_h */
