#!/usr/bin/env python3
#basic_functions.py


#Declare function to multiply two integers
def multiply(a, b):
    """Will return the product of two numbers, a and by
    """
    return(a*b)
 

#Declare a variable x to hold the product of the function
x = multiply(4, 2)
#Print the product of the function
print("The product of x is {}".format(x))


#Declare function that returns a string making the base value 
#of the string "you"
def hello_name(name="you"):
    """Will print the "name" string passed in the function
    """
    print("Hello, "+ name + "!")


#Pass two arguments to the function, a string and nothing to 
#demonstrate the two outcomes to the function
hello_name("Joseph")
hello_name()


#Declare function with myList as a variable
def less_than_ten(myList):
    """Will loop through the values of "myList" and return 
    only those less than 10
    """
    for n in myList:
        if n <= 10:
            """Print all the values less than 10 in myList"""
            print(n)


#Declare a list of values both less & greater than 10
number = [1,17,3,7,8,32]
#Pass the list of values as an argument to the less_than_ten function
less_than_ten(number)
