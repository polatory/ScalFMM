/* ceil example */
#include <stdio.h>      /* printf */
#include <math.h>       /* ceil */

int main ()
{
  printf ( "floor of 2.3 is %.1lf\n", floor (2.3) );
  printf ( "floor of 3.8 is %.1lf\n", floor (3.8) );
  printf ( "floor of -2.3 is %.1lf\n", floor (-2.3) );
  printf ( "floor of -3.8 is %.1lf\n\n", floor (-3.8) );
 
  printf ( "ceil of 2.3 is %.1f\n", ceil(2.3) );
  printf ( "ceil of 3.8 is %.1f\n", ceil(3.8) );
  printf ( "ceil of -2.3 is %.1f\n", ceil(-2.3) );
  printf ( "ceil of -3.8 is %.1f\n\n", ceil(-3.8) );

 printf ( "lround (2.3) = %ld\n", lround(2.3) );
  printf ( "lround (3.8) = %ld\n", lround(3.8) );
  printf ( "lround (-2.3) = %ld\n", lround(-2.3) );
  printf ( "lround (-3.8) = %ld\n", lround(-3.8) );

  return 0;


}
