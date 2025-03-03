# estimatelog
This is a parallel algorithm that computes the natural logarithm fucntion by numberical approximation. Numerical integration is used to compute the area under the curve, using the rectangle rule. 

### Build with
mpicc -Wall -o estimat_log estimate_log.c -lm

### Usage 
mpirun \[MPI options\] estimate_log <log_number> <num_intervals>
The user must supply and input number that the logarithm is to be determined of and a number of interval to be used in the approximation.

### Output and Features
*input_number*<tab>*computed_value*<tab>*absolute_error*<tab>*elapsed_time seconds*
The program outputs the input number followed by the computed logarithim value. Then the absolute error is outputed which is the absoulte value of the differnce between the log function and the computed value

### Defects/Shortcomings
-Originally all of the processors proccessed the command line, which isnt favorable but that was corrected with the use of MPI_BCAST.

Thank you to Professor Weiss for the assigment!
