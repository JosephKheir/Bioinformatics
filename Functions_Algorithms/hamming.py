#!/usr/bin/env python3
#hamming.py


#Import module sys
import sys


#Declare a function to find the hamming distance between two strings
def hamming(string1, string2):
    """Will return the hamming distance between two strings passed 
    into the function
    """
    dist = 0
    length = len(string1)
    for i in range(length):
        #Make sure both strings are equal in length
        if string1[i] != string2[i]:
            dist += 1
    return dist        


#Generate a code block to call the function from the CLI
if __name__=="__main__":
    #Make sure at least 2 arguments are passed
    arg_count = len(sys.argv) -1
    if arg_count < 2:
        #If less than 2 arguments are passed, raise exception
        raise Exception("This script requires at least 2 arguments")
    
    
    #Declare value of string 1 & 2 as argument 1 & 2 passed in the CLI
    string1 = sys.argv[1]
    string2 = sys.argv[2]


#Print hamming distance between two strings passed in CLI
print(hamming(string1, string2))
