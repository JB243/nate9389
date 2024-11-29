---
title: "(Template) Homework 4 Solution"
subtitle: "<h2><u>Bioinformatics 501 (Term Year)</u></h2>"
author: "<h3>Student Name</h3>"
date: "`r format(Sys.time(), '%B %Y')`"
tags: [DSPA, SOCR, MIDAS, Big Data, Predictive Analytics] 
output:
  html_document:
    theme: spacelab
    highlight: tango
    toc: true
    number_sections: true
    toc_depth: 2
    toc_float:
      collapsed: false
      smooth_scroll: true
    code_folding: show
    self_contained: yes
---


This is a template for my first HM problem submitted as RMD notebook
**Template_HW_4_Optimization_UMich_Bioinfo501.Rmd**.

# Module 4 (Optimization)

Review the R/RStudio installation and [R fundamentals](http://www.socr.umich.edu/people/dinov/courses/DSPA_notes/01_Foundation.html) as well as the [Optimization Chapter](http://www.socr.umich.edu/people/dinov/courses/DSPA_notes/21_FunctionOptimization.html) of the [DSPA Textbook](http://www.socr.umich.edu/people/dinov/courses).
 

## Problem 1
Find the extrema of this variant of the [Rosenbrock function](https://en.wikipedia.org/wiki/Rosenbrock_function), $(100(x^2−y)^2+(x−1)^2+100(z^2−w)^2+(z−1)^2)$, $N=4$, and compare your estimates against the Wolfram Alpha Optimizer.

## Solution 1
For the given structure of the Rosenbrock function, there are no local maxima as the function can diverge easily. Therefore, solving the problem requires focusing on finding both a global minimum and any local minima.

I optimized this function using the BFGS method as follows:

```{r}
# Define the Rosenbrock variant function
rosenbrock_variant <- function(vars) {
    x <- vars[1]
    y <- vars[2]
    z <- vars[3]
    w <- vars[4]
    f <- 100 * (x^2 - y)^2 + (x - 1)^2 + 100 * (z^2 - w)^2 + (z - 1)^2
    return(f)
}

# Initial guess for optimization
initial_guess <- c(x = -1.2, y = 1, z = -1.2, w = 1)

# Perform optimization using BFGS
result <- optim(
    par = initial_guess,
    fn = rosenbrock_variant,
    method = "BFGS",
    control = list(reltol = 1e-8)
)

# Display the optimization results
cat("Local Minimum Value: ", result$value, "\n")
cat("Point of Local Minimum: x =", result$par[1], 
    ", y =", result$par[2], 
    ", z =", result$par[3], 
    ", w =", result$par[4], "\n")
```

As a result, I got the following results:

```{bash}
Local Minimum Value:  7.668621e-08 
Point of Local Minimum: x = 0.9998085 , y = 0.9996164 , z = 0.9998002 , w = 0.9995997 
```

By doing the above multiple times, I could acquire only one local minimum. Also, I used `GenSA` package to find a global minimum as follows:

```{r}
library(GenSA)

# Define the Rosenbrock variant function
rosenbrock_variant <- function(vars) {
    x <- vars[1]
    y <- vars[2]
    z <- vars[3]
    w <- vars[4]
    f <- 100 * (x^2 - y)^2 + (x - 1)^2 + 100 * (z^2 - w)^2 + (z - 1)^2
    return(f)
}

result <- GenSA(
    par = c(x = -1.2, y = 1, z = -1.2, w = 1),
    fn = rosenbrock_variant,
    lower = c(-5, -5, -5, -5),
    upper = c(5, 5, 5, 5)
)

cat("Global Minimum Value:", result$value, "\n")
cat("Point of Global Minimum: x =", result$par[1], 
    ", y =", result$par[2], 
    ", z =", result$par[3], 
    ", w =", result$par[4], "\n")
```

As a result, I got the following results:

```{bash}
Global Minimum Value: 3.994835e-20 
Point of Global Minimum: x = 1 , y = 1 , z = 1 , w = 1
```

Thus, there is one global minimum and local minimum where x = 1, y = 1, z = 1, w = 1, and the value of 0.

When, referring to [WolframAlpha](https://www.wolframalpha.com/input?i=extrema+100%28x%5E2-y%29%5E2%2B%28x-1%29%5E2%2B100%28z%5E2-w%29%5E2%2B%28z-1%29%5E2), I found the same result.


## Problem 2
Find the minimum of the constrained [Mishra's Bird function](https://en.wikipedia.org/wiki/Test_functions_for_optimization#Test_functions_for_constrained_optimization), $f(x,y)=\sin(y)e^{[(1−\cos x)^2]}+\cos(x)e^{[(1−\sin y)^2]}+(x−y)^2$, subject to the constraint $(x+5)^2+(y+5)^2<25$.

Here we will be doing the WH4.2 solution and R code.

```{r}
# enter the R code here ....

```
  

## Prpoblem 3
Minimize the function $f(x,y,z)=−(x^3+5y−2^z)$, subject to 
$$\begin{cases}
x -\frac{y}{2}+z^2 \leq 50\\
\mod(x, 4) + \frac{y}{2} \leq 1.5
\end{cases} .$$

Check you solution against the Wolfram Alpha solution.

 
Here we will be doing the WH4.3 solution and R code.

```{r}
# enter the R code here ....

```
 

Your submission should include two files: RMD (source) and HTML (knitted report).
