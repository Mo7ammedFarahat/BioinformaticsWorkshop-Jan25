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

### Why Use R?
R is a powerful programming language and environment for statistical computing and graphics. It is widely used for data analysis, statistical modeling, and creating visually appealing charts. R is popular among statisticians, data scientists, and researchers for its versatility and large ecosystem of packages.

## Basic Syntax

### Arithmetic Operations
In R, arithmetic operations like addition, subtraction, multiplication, and division can be performed with simple syntax. These operations are fundamental in any programming language and are often the starting point for most calculations.

```R
# Basic arithmetic
2 + 3     # Addition
10 - 5    # Subtraction
4 * 6     # Multiplication
20 / 4    # Division
```

### Variables in R
Variables are used to store data values in programming. In R, you can assign values to variables using `<-` (preferred) or `=`. Variables hold data like numbers, strings, or logical values that you can manipulate.

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

## R Reserved Words

### Explanation of Reserved Words
Reserved words are predefined keywords in R that have special meanings and cannot be used as variable names. These words are used to perform specific actions like conditional checks, loops, or logical operations.

Examples include `if`, `else`, `repeat`, `while`, `TRUE`, `FALSE`, `NULL`, etc.

### Why Reserved Words are Important
Reserved words are integral to R programming because they help define control structures, logical conditions, and the flow of execution. Using these words as variable names would cause errors, as R wouldn't be able to distinguish between the reserved word and a variable.

## Types of Operators

### Arithmetic Operators
Arithmetic operators in R allow you to perform basic mathematical operations such as addition, subtraction, multiplication, and division.

```R
a <- 10
b <- 5
print(a + b)  # Output: 15
print(a - b)  # Output: 5
print(a * b)  # Output: 50
print(a / b)  # Output: 2
```

### Relational Operators
Relational operators are used to compare two values and return a logical result (`TRUE` or `FALSE`). These comparisons are essential for conditional logic and decision-making.

```R
a <- 10
b <- 5
print(a > b)   # Output: TRUE
print(a < b)   # Output: FALSE
print(a == b)  # Output: FALSE
print(a != b)  # Output: TRUE
```

### Logical Operators
Logical operators are used to manipulate logical values (`TRUE` and `FALSE`). These operators combine or invert logical conditions and are fundamental in conditional expressions.

```R
a <- TRUE
b <- FALSE
print(a & b)  # Output: FALSE (AND operation)
print(a | b)  # Output: TRUE (OR operation)
print(!a)     # Output: FALSE (NOT operation)
```

### Assignment Operators
Assignment operators are used to assign values to variables. The most common operator in R is `<-`, but `=` can also be used.

```R
a <- 10  # Assign value 10 to a
b <<- 20 # Global assignment
c = 30   # Assign value 30 to c
print(a)
print(b)
print(c)
```

### Operator Precedence
Operator precedence determines the order in which operations are performed in an expression. In mathematical expressions, operators like multiplication and division take precedence over addition and subtraction.

```R
result <- 10 + 5 * 2  # Multiplication happens first
print(result)         # Output: 20
```

## R Data Types

### Logical Data Type
In R, logical data types are used to represent truth values (`TRUE` or `FALSE`). They are essential for decision-making in conditional expressions.

```R
x <- TRUE       # Logical
y <- FALSE      # Logical
```

### Numeric Data Type
Numeric data types in R represent real numbers (decimals) and are commonly used for mathematical and statistical calculations.

```R
x <- 3.14       # Numeric
```

### Integer Data Type
In R, integers are whole numbers and are denoted with an `L` suffix. While numeric types can represent both integers and real numbers, integers are used when you specifically want to work with whole numbers.

```R
x <- 42L        # Integer
```

### Complex Data Type
Complex numbers consist of a real and an imaginary part. In R, complex numbers are represented as `a + bi`.

```R
x <- 1 + 2i     # Complex
```

### Character Data Type
Character data types store text data, which can be letters, words, or sentences. Strings are enclosed in double or single quotes.

```R
x <- "R Programming"  # Character
```

## R Data Object Types

### Vectors
A vector is the most basic data structure in R. It is a collection of elements of the same type. You can create a vector using the `c()` function.

```R
vec <- c(1, 2, 3, 4, 5)
print(vec)
```

### Matrices
A matrix is a two-dimensional array where each element has the same data type. You can create matrices using the `matrix()` function, specifying the number of rows and columns.

```R
matrix <- matrix(1:6, nrow=2, ncol=3)
print(matrix)
```

### Lists
A list is an ordered collection of objects that can be of different types. Lists are more flexible than vectors as they can store mixed data types.

```R
my_list <- list("Apple", 42, TRUE)
print(my_list)
```

### Arrays
An array is a multi-dimensional data structure that can hold values of the same data type. Unlike matrices, arrays can have more than two dimensions.

```R
arr <- array(1:8, dim = c(2, 2, 2))
print(arr)
```

### Factors
Factors are used to represent categorical data with fixed levels, such as "Low", "Medium", or "High". They are particularly useful for storing data that can be classified into distinct categories.

```R
factor_data <- factor(c("High", "Medium", "Low", "Medium"))
print(factor_data)
```

### Data Frames
A data frame is a two-dimensional table that can store data of different types in each column. Data frames are particularly useful for working with datasets.

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
Functions in R allow you to define reusable blocks of code. You can define your own functions using the `function()` keyword, which makes your code more modular and efficient.

```R
add <- function(a, b) {
  return(a + b)
}
print(add(5, 10))  # Output: 15
```

### Packages
R has a rich ecosystem of packages that extend its functionality. You can install and use these packages to access additional tools, functions, and datasets. For example, `ggplot2` is a popular package for creating visualizations.

```R
install.packages("ggplot2")
library(ggplot2)
```

## R File Manipulation

### CSV Files
CSV files (Comma Separated Values) are a common format for storing tabular data. You can read and write CSV files in R using built-in functions like `read.csv()` and `write.csv()`.

```R
data <- read.csv("file.csv")
print(data)
```

### Excel Files
R also supports reading Excel files using the `readxl` package. This allows you to work with spreadsheet data directly in R.

```R
install.packages("readxl")
library(readxl)
data <- read_excel("file.xlsx")
print(data)
```

## R Data Visualization

### Pie Charts
A pie chart is a circular visualization used to represent proportions of a whole. Each slice represents a category’s proportion.

```R
slices <- c(10, 20, 30)
labels <- c("A", "B", "C")
pie(slices, labels)
```

### Bar Charts
Bar charts represent data with rectangular bars, where the length of each bar corresponds to the value of the variable being plotted.

```R
barplot(c(10, 20, 30))
```

### Boxplots
Boxplots are used to visualize the distribution of numerical data by showing the median, quartiles, and potential outliers.

```R
boxplot(mpg ~ cyl, data=mtcars)
```

### Histograms
A histogram is used to represent the frequency distribution of a set of continuous or discrete data.

```R
hist(mtcars$mpg)
```

## Conclusion

This concludes the second day of our R course. Today, we've explored basic syntax, operators, data types, data structures, and visualizations. By the end of this day, you should be familiar with essential concepts in R and be ready to dive deeper into more advanced topics in the upcoming lessons.

