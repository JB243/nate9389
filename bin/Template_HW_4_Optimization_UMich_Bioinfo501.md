---
title: "Homework 4 Solution"
subtitle: "<h2><u>Bioinformatics 501 (Fall 2024)</u></h2>"
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


This is an RMD notebook **Template_HW_4_Optimization_UMich_Bioinfo501.Rmd**.

# Module 4 (Optimization)

Review the R/RStudio installation and [R fundamentals](http://www.socr.umich.edu/people/dinov/courses/DSPA_notes/01_Foundation.html) as well as the [Optimization Chapter](http://www.socr.umich.edu/people/dinov/courses/DSPA_notes/21_FunctionOptimization.html) of the [DSPA Textbook](http://www.socr.umich.edu/people/dinov/courses).
 

## Packages
-   R 4.4.0
-   alabama==2023.1.0
-   GenSA==1.1.14.1


## Problem 1
Find the extrema of this variant of the [Rosenbrock function](https://en.wikipedia.org/wiki/Rosenbrock_function), $100(x^2−y)^2+(x−1)^2+100(z^2−w)^2+(z−1)^2$, $N=4$, and compare your estimates against the Wolfram Alpha Optimizer.


## Solution 1
The structure of the given Rosenbrock function ensures that no local maxima exist, as the function tends to diverge. For example, if a particular variable approaches infinity while the other variables remain constant, the value of the function also tends toward infinity. Therefore, solving the problem requires focusing on finding both the global minimum and any local minima. However, it is evident that the function is differentiable at all points and its local minima do not occur at infinity for any variable. Therefore, the local minima of the function can be determined through an optimization process aimed at finding the minimum value(s).

Nonetheless, when updating $(x,y,z,w)$ in the direction that minimizes the function value, adjustments can always be made to satisfy $x^2=y$, $x=1$, $z^2=w$, and $z=1$. Hence, it is clear that there is only one minimum point, and that point is $x=y=z=w=1$, making the functional value of $0$.

I optimized the Rosenbrock function using various methods as follows:

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

# Optimization strategies
strats <- c("CG", "Nelder-Mead", "BFGS", "L-BFGS-B", "SANN")

# Set seed for reproducibility
set.seed(123)

# Perform optimization using all strategies
results <- lapply(strats, function(method) {
  control_params <- if (method == "L-BFGS-B") {
    list(factr = 1e7, pgtol = 1e-8)  
  } else {
    list(reltol = 1e-8) 
  }
  
  optim(
    par = initial_guess,
    fn = rosenbrock_variant,
    method = method,
    control = control_params
  )
})

# Display results
for (i in seq_along(strats)) {
  cat("Method:", strats[i], "\n")
  cat("  Local Minimum Value:", results[[i]]$value, "\n")
  cat("  Point of Local Minimum: x =", results[[i]]$par[1],
      ", y =", results[[i]]$par[2],
      ", z =", results[[i]]$par[3],
      ", w =", results[[i]]$par[4], "\n\n")
}
```

As a result, I got the following results:

```{bash}
Method: CG 
  Local Minimum Value: 0.01923853 
  Point of Local Minimum: x = 0.8996096 , y = 0.8084839 , z = 0.8996096 , w = 0.8084839 

Method: Nelder-Mead 
  Local Minimum Value: 0.03607516 
  Point of Local Minimum: x = 1.130008 , y = 1.276655 , z = 0.8615604 , w = 0.742373 

Method: BFGS 
  Local Minimum Value: 7.668621e-08 
  Point of Local Minimum: x = 0.9998085 , y = 0.9996164 , z = 0.9998002 , w = 0.9995997 

Method: L-BFGS-B 
  Local Minimum Value: 9.427752e-07 
  Point of Local Minimum: x = 1.000487 , y = 1.000976 , z = 0.9991604 , w = 0.9983197 

Method: SANN 
  Local Minimum Value: 0.01997143 
  Point of Local Minimum: x = 1.008519 , y = 1.017541 , z = 1.124792 , w = 1.258595  
```

Based on the results above, the different methods have varying efficacy based on the problem. With different seed values, I could acquire only one local minimum. Also, I used `GenSA` package to find a global minimum as follows:

```{r}
library(GenSA)

# Set seed for reproducibility
set.seed(123)

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

In this case, there is one global minimum and local minimum where $x=1, y=1, z=1, w=1$, with a function value of 0. Here, $3.994835\times 10^{-20}$ can be regarded as zero. At the beginning of the solution, it seems that various optimization methods yield similar results, differing only in how close they are to the solution. Ultimately, they all appear to converge to the same conclusion. 

Referring to the result from [WolframAlpha](https://www.wolframalpha.com/input?i=extrema+100%28x%5E2-y%29%5E2%2B%28x-1%29%5E2%2B100%28z%5E2-w%29%5E2%2B%28z-1%29%5E2), I confirmed the same conclusion.


## Problem 2
Find the minimum of the constrained [Mishra's Bird function](https://en.wikipedia.org/wiki/Test_functions_for_optimization#Test_functions_for_constrained_optimization), $f(x,y)=\sin(y)\exp{[(1−\cos x)^2]}+\cos(x)\exp{[(1−\sin y)^2]}+(x−y)^2$, subject to the constraint $(x+5)^2+(y+5)^2<25$.


## Solution 2

I optimized the Mishra's Bird function using the `alabama` package. 

```{r}
# Load required package
library(alabama)

# Define the objective function
f_obj <- function(xy) {
  x <- xy[1]
  y <- xy[2]
  sin(y) * exp((1 - cos(x))^2) +
  cos(x) * exp((1 - sin(y))^2) +
  (x - y)^2
}

# Define the constraint function
hin <- function(xy) {
  25 - (xy[1] + 5)^2 - (xy[2] + 5)^2
}

# Initial guess
x0 <- c(-3, -3)

# Set seed for reproducibility
set.seed(123)

# Perform optimization
result <- auglag(par = x0, fn = f_obj, hin = hin)

# Print results
cat("Global Minimum Value:", result$value, "\n")
cat("Point of Global Minimum: x =", result$par[1], 
    ", y =", result$par[2], "\n")
```
  
As a result, I got the following results:

```{bash}
Global Minimum Value: -106.7645 
Point of Global Minimum: x = -3.130308 , y = -1.582165
```

In this case, there is one global minimum and local minimum where $x=-3.130308, y=-1.582165$, with a function value of $-106.7645$. Here, the solution met the constraint given $(-3.130308+5)^2+(-1.582165+5)^2 = 15.17734 < 25$.

I also used `Rsolnp` package for verification.

```{r}
# Load the required package
library(Rsolnp)

# Define the objective function
f_obj <- function(xy) {
  x <- xy[1]
  y <- xy[2]
  sin(y) * exp((1 - cos(x))^2) +
    cos(x) * exp((1 - sin(y))^2) +
    (x - y)^2
}

# Define the inequality constraint function
ineq_fun <- function(xy) {
  # Constraint: (x + 5)^2 + (y + 5)^2 < 25
  # Transform to ineq_fun(xy) ≥ 0
  25 - (xy[1] + 5)^2 - (xy[2] + 5)^2
}

# Set the initial guess
x0 <- c(-3, -3)

# Set seed for reproducibility
set.seed(123)

# Perform optimization
result <- solnp(par = x0, 
                fun = f_obj, 
                ineqfun = ineq_fun, 
                ineqLB = 0, 
                ineqUB = Inf)

# Print the results
cat("Minimum Value:", result$values[length(result$values)], "\n")
cat("Point of Minimum: x =", result$pars[1], ", y =", result$pars[2], "\n")
```

As a result, I can acquire the following results:

```{bash}
Minimum Value: -106.7645 
Point of Minimum: x = -3.130247 , y = -1.582142
```

These results can be validated through visual inspection using the Plotly visualization.

Referring to the result from [WolframAlpha](https://www.wolframalpha.com/input?i=minimum+sin%28y%29+*+exp%28%281+-+cos%28x%29%29%5E2%29+%2B+++++cos%28x%29+*+exp%28%281+-+sin%28y%29%29%5E2%29+%2B++%28x+-+y%29%5E2%2C+subject+to+%28x%2B5%29%5E2+%2B+%28y%2B5%29%5E2+%3C+25), I confirmed the same outcome.


## Problem 3
Minimize the function $f(x,y,z)=−(x^3+5y−2^z)$, subject to $x -\frac{y}{2}+z^2 \leq 50$ and $\text{mod}(x, 4) + \frac{y}{2} \leq 1.5$.

Check your solution against the Wolfram Alpha solution.

 
## Solution 3

I referred to the source code presented in the lecture note. Note that $\text{mod}(a, n) = a - n \times \lfloor \frac{a}{n} \rfloor$ holds, e.g. $\text{mod}(4.31, 4) = 0.31$. 


```{r}
# Define the objective function
fun3 <- function (x){
  -(x[1]^3 + 5*x[2] - 2^x[3])
}

# Define the constrained objective function
fun3_constraints <- function(x) {
  if (x[1] - (x[2]/2)+x[3]^2 > 50) {NA}       # constraint 1
  else if ((x[1] %% 4) + (x[2]/2) > 1.5) {NA} # constraint 2
  else { fun3(x) } 
}

# Set initial parameters
x0 <- c(0.01, 2, -2)

# Initialize variables to store the best result
x_optimal_value <- Inf
x_optimal_point <- c(NA, NA)

# Set seed for reproducibility
set.seed(123)

# Run 50 iterations 
for (i in 1:50) {
    x_optimal <- optim(x0, fun3_constraints)
    if (x_optimal$value < x_optimal_value) {
        x_optimal_value <- x_optimal$value
        x_optimal_point <- x_optimal$par
    }
}

# Display the results
cat("The minimum value of the function is:", x_optimal_value, "\n")
cat("This occurs at x =", x_optimal_point[1], ", y =", x_optimal_point[2], ", z =", x_optimal_point[3], "\n")
```

As a result, I got the following results:

```{bash}
The minimum value of the function is: -14.76122 
This occurs at x = 7.225925e-09 , y = 3 , z = -2.066219 
```

Note that the $x -\frac{y}{2}+z^2 \leq 50$ holds, but $\text{mod}(x, 4) + \frac{y}{2} \leq 1.5$ doesn't hold. However, due to the nature of numerical algorithm, the obtained solution is an approximate value, and if $x$ is positive, $y$ would take a value slightly less than 3, thus satisfying the given condition.

I also used simulated annealing method, which is a slower stochastic global optimization optimizer that works well with difficult functions, e.g., non-differentiable, non-convex. 

```{r}
# Define the objective function
fun3 <- function (x){
  # Objective function to minimize
  return(-(x[1]^3 + 5*x[2] - 2^x[3]))
}

# Define the constrained objective function
fun3_constraints <- function(x) {
  # Constraint 1: x - y/2 + z^2 <= 50
  if (x[1] - (x[2]/2) + x[3]^2 > 50) {
    return(NA)  # Penalize if constraint is violated
  }
  # Constraint 2: mod(x, 4) + y/2 <= 1.5
  if ((x[1] %% 4) + (x[2]/2) > 1.5) {
    return(NA)  # Penalize if constraint is violated
  }
  # If all constraints are satisfied, return the objective function value
  return(fun3(x))
}

# Set initial parameters
x0 <- c(0.01, 2, -2)

# Initialize variables to store the best result
best_value <- Inf
best_par <- NULL

# Set seed for reproducibility
set.seed(123)

# Run multiple iterations to find the global minimum
for (i in 1:50) {
  result <- optim(
    par = x0,
    fn = fun3_constraints,
    method = "SANN",
    control = list(maxit = 10000)
  )
  # Update the best result if current result is better
  if (!is.na(result$value) && result$value < best_value) {
    best_value <- result$value
    best_par <- result$par
  }
}

# Display the results
cat("The minimum value of the function is:", best_value, "\n")
cat("This occurs at x =", best_par[1], ", y =", best_par[2], ", z =", best_par[3], "\n")
```

As a result I met the following result: 

```{bash}
The minimum value of the function is: -108839.3 
This occurs at x = 47.74851 , y = -4.49735 , z = 0.02309909 
```

Note that the $x -\frac{y}{2}+z^2 \leq 50$ and $\text{mod}(x, 4) + \frac{y}{2} \leq 1.5$ both hold. 

The difference between two is the optimization method used; that is, the first method uses a default Nelder-Mead method while the second method uses SANN (simulated annealing) method. However, I have also observed cases where SANN converges to ~ -14.8 instead of ~ -109000. The geometry of the objective function includes stable local minima and multiple divergent points, which means that depending not only on the optimization algorithm but also on the initial starting point and search path, the solution may converge to a local minimum or become trapped by the boundaries of the geometry. This implies that optimization algorithms do not guarantee exploration of the entire solution space. SANN (Simulated Annealing) tends to focus more on exploring the global minimum, making it more likely to produce an improved solution.

I created a visual representation of a grid of points that adhere to the constraints, with each point colored based on its corresponding function value. While examining the plot, I noticed areas with even lower function values. For exameple, the point $x=49.5, y=0, z=0$ meets the specified constratints and yields a function value of $-121286.4$.

Referring to the result from [WolframAlpha](https://www.wolframalpha.com/input?i=minimize+-%28x%5E3+%2B+5*y+-2%5Ez%29++subject+to+x+-0.5*y%2Bz%5E2+%3C%3D+50+AND+mod%28x%2C+4%29+%2B+y%2F2+%3C%3D+1.5), I observed that the result aligns with the initial outcome mentioned above. However, considering other lower function values, it suggests that both WolframAlpha and the optimization algorithm might have inaccuracies.
