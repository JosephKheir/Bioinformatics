#!/usr/bin/env python3
#hello_you.py


#Import the module sys
import sys


#Declare the function with a string variable "name"
def hello_name(name):
    """Will return the value of the string passed 
    to the function
    """
    return("Hello, {name}!")


#Generate a code block to call the function from the CLI
if __name__ == '__main__':
    #Make sure at least one argument is passed 
    arg_count = len(sys.argv) - 1
    if arg_count !=1:
        #If only one argument is passed, raise and exception
        #to exit the script
        raise Exception("Hello, you!")


    #Declare a variable to house the value of the first argument
    #raised in the CLI and pass to the function
    x = str(sys.argv[1])


    #Print the value of the function
    print("Hello, {}!".format(x))

