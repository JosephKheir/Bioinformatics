#!/usr/bin/env python
# add_numbers.py

# Initialize two totals variables
total = 0
total2 = 0

# Iterate over a range of values
for num in range(1,10):
    # Add num to total and assign to total
    total = total + num
    # Or use this shortcut to add and assign to total in one step
    total2 += num

print(total)
print(total2)

