/* Error management system
 */

#include<stdio.h>
#include"locerror.h"

void
locerror(char *function, char *message)
{
	fprintf(stderr, "Fatal error : function %s, %s\n", function, message);
	exit(1);
}
