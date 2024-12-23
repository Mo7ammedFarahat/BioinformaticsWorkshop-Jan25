* [Day2: Introduction to R programming language](#day2-introduction-to-r-programming-language)
    - [R & R Studio Installing](#r--r-studio-installing)
    - [Basic Syntax](#basic-syntax)
    - [R Reserved Words](#r-reserved-words)
    - [Variables](#variables)
    - [Types of Operators](#types-of-operators)  
        - [Arithmetic Operators](#arithmetic-operators)  
        - [Relational Operators](#relational-operators)  
        - [Logical Operators](#logical-operators)  
        - [Assignment Operators](#assignment-operators)
    - [Operator Precedence](#operator-precedence)
    - [R Data Types](#r-data-types)  
        - [Logical](#logical)  
        - [Numeric](#numeric)  
        - [Integer](#integer)  
        - [Complex](#complex)  
        - [Character](#character)
    - [R Data Object Types](#r-data-object-types)  
        - [Vectors](#vectors)  
        - [Matrices](#matrices)  
        - [Lists](#lists)  
        - [Arrays](#arrays)  
        - [Factors](#factors)  
        - [Data Frames](#data-frames)
    - [R Code Blocks](#r-code-blocks) 
        - [R - Functions](#r---functions)
        - [R - Packages](#r---packages)
    - [R File manipulation](#r-file-manipulation) 
        - [R - CSV Files](#r---csv-files)
        - [R - Excel File](#r---excel-file)
    - [R Data Visualization](#r-data-visualization) 
        - [R - Pie Charts](#r---pie-charts)
        - [R - Bar Charts](#r---bar-charts)
        - [R - Boxplots](#r---boxplots)
        - [R - Histograms](#r---histograms)
        - [R - Line Graphs](#r---line-graphs)
        - [R - Scatterplots](#r---scatterplots)
        - [R - ggplot](#r---ggplot)
    - [R Basic Statistics](#r-basic-statistics) 
        - [R - Mean, Median and Mode](#r---mean-median-and-mode)
        - [R - Normal Distribution](#r---normal-distribution)
        - [R - Binomial Distribution](#r---binomial-distribution)

# Day 2: Introduction to R Programming Language

Welcome to Day 2 of the R course! This section provides an in-depth introduction to the R programming language, covering installation, syntax, data manipulation, visualization, and basic statistics.

## R & R Studio Installing

To start using R, you need to install both R and RStudio.  
1. **Installing R**: Visit [The Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/), choose your operating system, and follow the installation instructions.  
2. **Installing RStudio**: Visit [RStudio](https://www.rstudio.com/), download the IDE, and install it.

## Basic Syntax

R uses a simple syntax for programming. Here's an example:

```R
# Basic arithmetic
2 + 3     # Addition
10 - 5    # Subtraction
4 * 6     # Multiplication
20 / 4    # Division

# Assigning values to variables
x <- 10
y <- 5
z <- x + y
print(z)  # Output: 15
```

## R Reserved Words

Reserved words are predefined keywords in R that you cannot use as variable names. Examples include `if`, `else`, `repeat`, `while`, `function`, `TRUE`, `FALSE`, `NULL`, etc.

## Variables

Variables in R store data values. You can assign values to variables using `<-` or `=`.

```R
# Assigning values
name <- "R Programming"
age <- 25
isStudent <- TRUE

# Displaying variable values
print(name)       # Output: "R Programming"
print(age)        # Output: 25
print(isStudent)  # Output: TRUE
```

## Types of Operators

R has several types of operators to perform different operations.

### Arithmetic Operators

```R
a <- 10
b <- 5
print(a + b)  # Output: 15
print(a - b)  # Output: 5
print(a * b)  # Output: 50
print(a / b)  # Output: 2
```

### Relational Operators

```R
a <- 10
b <- 5
print(a > b)   # Output: TRUE
print(a < b)   # Output: FALSE
print(a == b)  # Output: FALSE
print(a != b)  # Output: TRUE
```

### Logical Operators

```R
a <- TRUE
b <- FALSE
print(a & b)  # Output: FALSE
print(a | b)  # Output: TRUE
print(!a)     # Output: FALSE
```

### Assignment Operators

```R
a <- 10  # Assign value 10 to a
b <<- 20 # Global assignment
c = 30   # Assign value 30 to c
print(a)
print(b)
print(c)
```

## Operator Precedence

Operator precedence determines the order of operations in an expression.

```R
result <- 10 + 5 * 2  # Multiplication happens first
print(result)         # Output: 20
```

## R Data Types

R supports various data types, such as logical, numeric, integer, complex, and character.

```R
x <- TRUE       # Logical
y <- 3.14       # Numeric
z <- 42L        # Integer
w <- 1+2i       # Complex
name <- "R"     # Character
```

## R Data Object Types

### Vectors

```R
vec <- c(1, 2, 3, 4, 5)
print(vec)
```

### Matrices

```R
matrix <- matrix(1:6, nrow=2, ncol=3)
print(matrix)
```

### Lists

```R
my_list <- list("Apple", 42, TRUE)
print(my_list)
```

### Arrays

```R
arr <- array(1:8, dim = c(2, 2, 2))
print(arr)
```

### Factors

```R
factor_data <- factor(c("High", "Medium", "Low", "Medium"))
print(factor_data)
```

### Data Frames

```R
df <- data.frame(
  Name = c("Alice", "Bob"),
  Age = c(25, 30),
  Gender = c("F", "M")
)
print(df)
```

## R Code Blocks

### Functions

```R
add <- function(a, b) {
  return(a + b)
}
print(add(5, 10))
```

### Packages

```R
install.packages("ggplot2")
library(ggplot2)
```

## R File Manipulation

### CSV Files

```R
data <- read.csv("file.csv")
print(data)
```

### Excel Files

```R
install.packages("readxl")
library(readxl)
data <- read_excel("file.xlsx")
print(data)
```

## R Data Visualization

### Pie Charts

```R
slices <- c(10, 20, 30)
labels <- c("A", "B", "C")
pie(slices, labels)
```

### Bar Charts

```R
barplot(c(10, 20, 30))
```

### Boxplots

```R
boxplot(mpg ~ cyl, data=mtcars)
```

### Histograms

```R
hist(mtcars$mpg)
```

### Line Graphs

```R
plot(1:10, type="l")
```

### Scatterplots

```R
plot(mtcars$mpg, mtcars$wt)
```

### ggplot

```R
library(ggplot2)
ggplot(mtcars, aes(x=mpg, y=wt)) + geom_point()
```

## R Basic Statistics

### Mean, Median, and Mode

```R
data <- c(1, 2, 3, 4, 5)
print(mean(data))
print(median(data))
```

### Normal Distribution

```R
x <- rnorm(100)
hist(x)
```

### Binomial Distribution

```R
x <- rbinom(100, size=10, prob=0.5)
hist(x)
```
