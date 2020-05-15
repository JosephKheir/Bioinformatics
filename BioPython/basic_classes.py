#!/usr/bin/env python3
#basic_classes.py


#Declare a class named Circle with attributes color and radius
class Circle():
    def __init__(self, color, radius):
        """ A function of the class circle with attributes color
        and radius, which are used in the proceeding functions
        to calculate and return diameter, circumfrence and color
        """
        self.color = color
        self.radius = radius
    
    
    #Declare function for Circle diameter
    def diameter(self):
        """ Calculates the diameter of a circle
        """
        return (self.radius * 2)
    
    
    #Declare function for Circle circumfrence
    def circumfrence(self):
        """ Calculates the circumfrence of a circle
        """
        return 2 * (3.14 * self.radius)

    
    #Declare function to check the color of Circle
    def isRed(self):
        """ Returns "True" if the circle color is red and 
        "False" if the circle color is not red
        """
        if(self.color == "red"):
            return True
        else:
            return False


#Declare a class named GraduateStudent
class GraduateStudent():
    def __init__(self, first_name, last_name, year, major):
        """ A function with attributes first name, last name, year
        and major, returns student first & last name, matriculation
        year and graduate major
        """
        self.first_name = first_name
        self.last_name = last_name
        self.year = year
        self.major = major

    
    def year_matriculated(self):
        """ Calculates year matriculated based off current year 
        of 2020
        """
        return (2020 - self.year)


#Generate variable to test functions of class Circle and GraduateStudent
line_separator = "----------"
circle1 = Circle("red", 4)
circle2 = Circle("pink", 17)
student1 = GraduateStudent("Joseph", "Kheir", 1, "BioInformatics")
student2 = GraduateStudent("Chauncey", "Lacey", 4, "Business Administration")


#Print results
print(f"circle1 color is red:", circle1.isRed())
print(f"circle1 diameter is:", circle1.diameter())
print(f"circle1 circumfrence is:", circle1.circumfrence())
print(line_separator)
print(f"circle2 color is red:", circle2.isRed())
print(f"circle2 diameter is:", circle2.diameter())
print(f"circle2 circumfrence is:", circle2.circumfrence())
print(line_separator)
print(student1.first_name, student1.last_name, student1.year_matriculated(), student1.major)
print(line_separator)
print(student2.first_name, student2.last_name, student2.year_matriculated(), student2.major)
