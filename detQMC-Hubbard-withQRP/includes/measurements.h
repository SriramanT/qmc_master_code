//
//  measurements.h
//  
//
//  Created by Francisco Brito on 10/07/2018.
//

#ifndef measurements_h
#define measurements_h

template<int L, int N>
class Measurements
{
    Eigen::Matrix<double, L, N> HSfield;
public:
    void genHsMatrix();
    Eigen::Matrix<double, L, N> matrix();
    double get(int x, int y);
    void flip(int l, int i);
};

template<int L, int N>
void Configuration<L, N>::genHsMatrix()
{
    //  Generate the HS field matrix
    int l;
    int i;
    HSfield = Eigen::Matrix<double, L, N>::Random(L,N);
    for (l = 0; l < L; l++)
    {
        for (i = 0; i < N; i++)
        {
            if ( HSfield(l, i) < 0 )
            {
                HSfield(l, i) = -1;
            }
            else
            {
                HSfield(l, i) = 1;
            }
        }
    }
}

#endif /* measurements_h */
