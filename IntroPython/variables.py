#!/usr/bin/env python
# variables.py
# A variable containing a string
name = 'Joseph Kheir'

# A variable containing an integer
age = 28

# A list of strings
names = ["Joseph Kheir", "Chauncey Lacey", "Veronica Wells"]

# A list of numbers
numbers = [15.5,  25.6, 50, 500]

# A dictionary of names and ages
ages = {'Joseph': 28, 'Chauncey': 31, 'Veronica': 56}

# Print the values of the variables
print(name)
print(age)
print(names)
print(numbers)
print(ages)

# Find out what type Python assigned to the variables
print(type(name))
print(type(age))
print(type(names))
print(type(numbers))
print(type(ages))

# Iterate over names in a for loop
for name in names:
	print(name)
	print(type(name))
