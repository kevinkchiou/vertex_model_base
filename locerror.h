/* Error management system
 */

#ifndef SRC_LOCERROR_LOCERROR_H
#define SRC_LOCERROR_LOCERROR_H 1

/* Exit the program with a message displaying the function where the error
 * occured with a message.
 */
void locerror(char *function, char *message);


void exit(int code);

#endif /* SRC_LOCERROR_LOCERROR_H */
