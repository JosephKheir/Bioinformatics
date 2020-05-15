#!/usr/bin/env python3
#fibonacci.py


#Import module sys
import sys


#Declare a function that calculates the population of C. elegans at a
#reproduction rate k on day n
def population(n, k):
    """Will return Population of C. elegans on day n
    """
    #Declare first fibonacci numbers as 0 and 1
    fib1, fib2 = 1, 1
    fib = 1
    for i in range(2, n):
        fib = fib1 + fib2 * k
        fib2, fib1 = fib1, fib
    return fib


#Generate a code block to call the function from the CLI
if __name__=="__main__":
    #Make sure at least 2 arguments are passed
    arg_count = len(sys.argv) -1
    if arg_count < 2:
        #If less than two arguments raise exception
        raise Exception("This script requires at least 2 arguments")
    
    
    #Declare input variables passed from the CLI
    n = int(sys.argv[1])
    k = int(sys.argv[2])


#Print population of C. elegans
print("The population after {} days is:".format(n))
print(population(n,k))
