#include <iostream>
#include "f1.hpp"
#include "Jacobian.hpp"
int main()
{
    otherdata junk;
    mat x={{0,2.6,3.5}};
    mat (*p)(mat, otherdata)=f1 ;
    std::cout<<Jacobian(p,x,junk);
    return 0;
}
