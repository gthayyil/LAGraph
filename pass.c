
#include <stdio.h>

void f (float *x)
{
    *x = *x + 1 ;
    printf ("x %g\n", *x) ;
}

int main (void)
{
    float x = 3.7 ;
    f (&x) ;
    printf ("x %g\n", x) ;
}

