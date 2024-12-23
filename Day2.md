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

A vector is a sequence of data elements of the same type.

```R
vec <- c(1, 2, 3, 4, 5)
print(vec)
```

### Matrices

A matrix is a two-dimensional array, where each element has the same data type.

```R
matrix <- matrix(1:6, nrow=2, ncol=3)
print(matrix)
```

### Lists

A list can hold elements of different types, such as vectors, numbers, or even other lists.

```R
my_list <- list("Apple", 42, TRUE)
print(my_list)
```

### Arrays

An array is similar to a matrix but can have more than two dimensions.

```R
arr <- array(1:8, dim = c(2, 2, 2))
print(arr)
```

### Factors

Factors are used to represent categorical data and store the levels of a variable.

```R
factor_data <- factor(c("High", "Medium", "Low", "Medium"))
print(factor_data)
```

### Data Frames

A data frame is a table or a two-dimensional array-like structure in which each column can contain different types of data.

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

Functions are used to group a set of instructions that can be called repeatedly.

```R
add <- function(a, b) {
  return(a + b)
}
print(add(5, 10))
```
Example 1: Working with DNA Sequences (Vectors)

In bioinformatics, sequences (like DNA or protein sequences) are often represented as character vectors. Here’s how we can use basic R functions to analyze DNA sequences.  

Example 1: Counting Nucleotides in a DNA Sequence

Let’s say we have a DNA sequence stored as a vector of characters:

```R
# DNA sequence as a vector of nucleotides
dna_sequence <- c("A", "T", "G", "C", "A", "A", "T", "G", "C", "G", "T", "A", "A")

# Counting the occurrences of each nucleotide
nucleotide_counts <- table(dna_sequence)
print(nucleotide_counts)
```
Explanation:
    `table(dna_sequence)` counts the number of occurrences of each nucleotide (A, T, G, C) in the sequence. This is useful in bioinformatics to determine the composition of a DNA sequence.  
<details>
  <summary>Output</summary>

```R
dna_sequence
A B C G
5 4 3 3
```
</details>
Example 2: Calculating GC Content of a DNA Sequence

The GC content is important in bioinformatics to understand the stability of a DNA strand. We can calculate the percentage of G and C bases in a DNA sequence.
```R
# Function to calculate GC content
calculate_gc_content <- function(dna_seq) {
  gc_count <- sum(dna_seq == "G" | dna_seq == "C")
  gc_content <- (gc_count / length(dna_seq)) * 100
  return(gc_content)
}

# DNA sequence
dna_seq <- c("A", "T", "G", "C", "A", "A", "T", "G", "C", "G", "T", "A", "A")

# Calculate GC content
gc_content <- calculate_gc_content(dna_seq)
print(paste("GC Content: ", gc_content, "%"))
```
<details>
  <summary>Output</summary>

  GC Content:  38.46 %

</details>

Example 3: Sequence Motif Search in a DNA Sequence
```R
# Define a DNA sequence and motif
sequence <- DNAString("ATCGATCGATCGATCG")
motif <- "ATCG"
```
<details>
  <summary>Answer</summary>

#### Search for the motif in the sequence
 `motif_positions <- vmatchPattern(motif, sequence)`

#### Print the positions of the motif
`motif_positions`  
Output:
`[1] 1 5 9 13`

</details>



### Packages

Packages in R extend the functionality of the base language.

```R
install.packages("ggplot2")
library(ggplot2)
```

## R File Manipulation

### CSV Files

You can read and write CSV files using `read.csv()` and `write.csv()`.

```R
data <- read.csv("file.csv")
print(data)
```

### Excel Files

To read Excel files, you can use the `readxl` package.

```R
install.packages("readxl")
library(readxl)
data <- read_excel("file.xlsx")
print(data)
```

## R Data Visualization

### Pie Charts

Pie charts are used to show the proportions of a whole.

```R
slices <- c(10, 20, 30)
labels <- c("A", "B", "C")
pie(slices, labels)
```

### Bar Charts

Bar charts represent categorical data with rectangular bars.

```R
barplot(c(10, 20, 30))
```

### Boxplots

Boxplots visualize the distribution of a dataset, highlighting the median and outliers.

```R
boxplot(mpg ~ cyl, data=mtcars)
```

### Histograms

Histograms represent the distribution of numerical data by dividing it into intervals.

```R
hist(mtcars$mpg)
```

### Line Graphs

Line graphs are used to represent data points in a time sequence, or any other continuous data. They are particularly useful for visualizing trends over time or ordered categories.

```R
# Example of a simple line graph
x <- c(1, 2, 3, 4, 5)
y <- c(2, 4, 6, 8, 10)
plot(x, y, type = "o", col = "blue", xlab = "X-axis", ylab = "Y-axis", main = "Simple Line Graph")
```

### Scatterplots

A scatterplot is a diagram that represents the relationship between two variables by plotting data points. It helps to identify patterns, trends, or correlations between variables.

```R
# Example of a scatterplot
x <- c(1, 2, 3, 4, 5)
y <- c(2, 4, 6, 8, 10)
plot(x, y, col = "red", pch = 19, xlab = "X-axis", ylab = "Y-axis", main = "Simple Scatterplot")
```

### ggplot

`ggplot2` is a powerful R package used for creating complex plots with simple code. It follows a grammar of graphics, allowing users to layer components of a plot.

```R
# Example using ggplot2
library(ggplot2)
data(mtcars)
ggplot(mtcars, aes(x = mpg, y = hp)) + 
  geom_point() + 
  ggtitle("Scatterplot of MPG vs Horsepower") +
  xlab("Miles Per Gallon (MPG)") + 
  ylab("Horsepower")
```

## R Basic Statistics

### Mean, Median, and Mode

- **Mean**: The average of a set of numbers.
- **Median**: The middle value when the numbers are sorted in ascending order.
- **Mode**: The most frequently occurring value in a dataset.

```R
# Example data
data <- c(1, 2, 2, 3, 4, 5, 6)

# Mean
mean_value <- mean(data)
print(mean_value)  # Output: 3.57

# Median
median_value <- median(data)
print(median_value)  # Output: 3

# Mode (custom function)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
mode_value <- getmode(data)
print(mode_value)  # Output: 2
```
`set.seed()` is used to set the starting point for generating random numbers. This ensures that the sequence of random numbers generated is reproducible. For example, if you use the same seed number, you will get the same random values every time you run the code. This is useful for reproducibility in experiments or simulations.

#### Finding Most Frequent Gene Expression (Mode Calculation)

In gene expression analysis, the mode can be useful to find the most frequently expressed value across different conditions or samples. Let’s calculate the mode using randomly generated data.

```R
# Set a seed to ensure reproducibility of random values
set.seed(789)

# Simulate random gene expression data (random values)
gene_expression <- sample(10:20, 100, replace = TRUE)
```
<details>
  <summary>Code</summary>

```R
# Function to calculate mode
calculate_mode <- function(x) {
  uniq_x <- unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}

# Calculate and print the mode
mode_value <- calculate_mode(gene_expression)
cat("Most Frequent Gene Expression (Mode):", mode_value, "\n")
```
#### Output:
`Most Frequent Gene Expression (Mode): 14`


### Normal Distribution

A **normal distribution** is a bell-shaped curve that represents the distribution of a dataset, where most of the data points are clustered around the mean. In R, you can simulate a normal distribution using the `rnorm()` function.

```R
# Generating a normal distribution
set.seed(123)
normal_data <- rnorm(1000, mean = 0, sd = 1)

# Visualizing the distribution
hist(normal_data, breaks = 30, col = "skyblue", main = "Normal Distribution", xlab = "Values")
```

### Binomial Distribution

A **binomial distribution** represents the number of successes in a fixed number of independent Bernoulli trials, each with the same probability of success.

```R
# Generating a binomial distribution
set.seed(123)
binomial_data <- rbinom(1000, size = 10, prob = 0.5)

# Visualizing the distribution
hist(binomial_data, breaks = 30, col = "green", main = "Binomial Distribution", xlab = "Number of Successes")
```

## Conclusion

In this continuation of Day 2, you’ve learned about various data visualization techniques such as line graphs, scatterplots, and using `ggplot2`. We also covered basic statistics concepts like mean, median, mode, normal distribution, and binomial distribution. These fundamental concepts will help you analyze and visualize data efficiently in R.
