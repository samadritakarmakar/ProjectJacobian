#ifndef JACOBIAN_HPP
#define JACOBIAN_HPP
#include <armadillo>

using namespace arma;
template<class OtherData>
mat Jacobian(mat (*f)(mat , OtherData), mat x1, OtherData data) //pointer of the function; and the column of vectors
{
    int i;
    vec x=vectorise(x1); //Ensures a column vector
    vec x_h=sqrt(eps(x))%x; //find x_h
    mat a;
    a<<1.0<<endr;
    double eps1=sqrt(eps(a).at(0,0));
    x_h.replace(0.0,eps1); //find elements where x_h=0 and replace them.
    //x_h.replace(0.0,.000000002); //fix if the above lines do not give good results.
    mat xPlus_h=kron(ones(x.n_rows,1),x.t())+diagmat(x_h); //Creates a matrix [x+x_h, y, z; x, y+y_h, z; x, y, z+z_h]
    const vec fconst= vectorise((*f)(x, data)); //values for f(x)
    vec fx1=vectorise((*f)(xPlus_h.row(0), data))-fconst; //for the 1st row; (f(x+x_h) - f(x))
    mat df_dx=zeros(fx1.n_rows, x.n_rows);
    for (i=0; i<x.n_rows-1; i++)
    {
        df_dx.col(i)=fx1/x_h(i); //(f(x+x_h)(i) - f(x))/x_h(i)
        fx1=vectorise((*f)(xPlus_h.row(i+1), data))-fconst;
    }
    df_dx.col(i)=fx1/x_h(i); //for the last column

    return df_dx;
}

#endif // JACOBIAN_HPP
