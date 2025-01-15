## Module.1: Bash

### Task 1
Write a Bash script that uses a `for` loop to list all `.fasta` files in the current directory and display the first three lines of each `.fasta` file, and also store the output in a new file called results.log.

<details>
<summary>Answer</summary>  

**Script:**
```bash
#!bin/bash

for file in *.fasta
do
  echo "these are the first 3 lines of $file:"
  head -n 3 $file
done > results.log
```
</details>

### Task 2 
Write a Bash script that checks if a file named `genome.fasta` exists in the current directory. If it exists, print "File found" else, print "File not found"

<details>
<summary>Answer</summary>

```bash
#!bin/bash

if [ -f "genome.fasta" ]; then
    echo "File found!"
else
    echo "File not found."
fi
```

</details>

---


## Module 2: R scripting

### Task 1
 
Create a data frame with for a gene expression matrix. which contains three genes `GeneA`, `GeneB`, `GeneC` in the first colmn, and theexpression values in the second colmn and print this data frame.

<details>
<summary>Answer</summary>

**Answer**:
```R
gene_data <- data.frame(
  Gene = c("GeneA", "GeneB", "GeneC"),
  Expression = c(10, 20, 15)
)

print(gene_data)
```

**Output**:
```
   Gene Expression
1 GeneA         10
2 GeneB         20
3 GeneC         15
```
</details>
---


### Task 2
Write a simple function that takes a numeric vector as input and returns the sum of its elements

<details>
<summary>Answer</summary>


**Answer**:
```R
sum_vector <- function(x) {
  return(sum(x))
}

Myvector <- c(1, 2, 3, 4)
result <- sum_vector(Myvector)

print(result)
```

**Output**:
```
10
```

</details>
