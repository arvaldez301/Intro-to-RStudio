# Intro-to-R
Introduction to Beginner R Programing. Here you will learn everything from basic R commands! Below you will find some steps to help get you comfortable with R Studion
## R Studio
To begin workign with R studio there are two programs that you must download. One being the core R software and the other being R studio where you wil conduct all of your R coding.

http://cran.r-project.org/

http://rstudio.com

R Studio also allows you to create Python scripts and run Terminal.
## Basics
### Arithmetic with R
Within R, you are able to preform basic mathematical functions.
```
# Calculate 3 + 4
3 + 4

# Calculate 6 + 12
6 + 12

# An addition
5 + 5 

# A subtraction
5 - 5 

# A multiplication
3 * 5

 # A division
(5 + 5) / 2 

# Exponentiation
2^5

# Modulo
28 %% 6
```
### Variable Assignment
Variable is a basic concept in statistical programing. A variable allows you to store a value (eg. 4) or an object (eg. a function description) in R. This value can later be used to easily access the value or the objet that is stored within this variable
```
my_var <_ 4
my_var
```
You can also combine multiple variables into a joined variable
```
# Assign a value to the variables my_apples and my_oranges
my_apples <- 5
my_oranges <- 6

# Add these two variables together
my_apples + my_oranges

# Create the variable my_fruit
my_fruit <- my_apples + my_oranges
```
### Basic Data Types in R
R is capable of workign with numerous data types some of the most basic ones to get started with are
  - Demical values (```4.5```) are called numerics
  - Whole numbers (```4```) are called integers, but are also numerics
  - Boolean values (```TRUE``` or ```FALSE```) are called logical
  - Text (or string) values are called characters

Quotation marks around text indications that is a string (```"some text"```)

If you wish to double check what data type you are working with you can call the ```class()``` function
```
# Declare variables of different types
my_numeric <- 42
my_character <- "universe"
my_logical <- FALSE 

# Check class of my_numeric
class(my_numeric)

# Check class of my_character
class(my_character)

# Check class of my_logical
class(my_logical)
```
## Vectors
### Creating a Vector
### Naming a Vector
