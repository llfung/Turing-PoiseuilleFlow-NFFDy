# Welcome to Julia!

# Julia is a high-level, high-performance programming language for technical computing.
# It was designed to be easy to use and at the same time provide the speed and power of traditional compiled languages.

# In this introductory material, we will cover some basic concepts and syntax of Julia.

## Running Julia Code in VS Code and the REPL
# Much like MATLAB, Julia can be run in an interactive environment called the REPL (Read-Eval-Print Loop), which is equivalent to the MATLAB command window.

# You can also run Julia code within the .jl file.
# To run the code in this file, click the Play button in the top right corner.
# To run a specific line of code, Press Ctrl + Enter.
# To run a block of code (separated by section header with ## in the beginning), Press Shift + Enter.

# There is a workspace in the VS Code that shows the variables and their values, similar to MATLAB.

## Variables and Data Types
# Unlike MATLAB and more like C, Julia require you to be explicit about the type of the variable
# - Integers: Int8, Int16, Int32, Int64 (default), Int128
# - Floating-point numbers: Float16, Float32, Float64 (default)

# Here's an example of declaring a variable and assigning a value to it:
x = 10 # x is an int64
y = 10.0 # y is a float64

# Variable names
# JULIA support unicode characters in variable names, which means you can use Greek letters and other symbols in your variable names.
# To type them, you can use LaTeX-like syntax. For example, to type the Greek letter α, you can type \alpha and then press the Tab key.

## Control Flow

# Julia provides various control flow statements, such as if-else, for loops, and while loops.

# Here's an example of an if-else statement:
if x > 5
    println("x is greater than 5")
else
    println("x is less than or equal to 5")
end

# Here's an example of a for loop:
for i in 1:5
    println(i)
end

# Vectors and Matrices
# JULIA's treatment to vectors and matrices is slightly different from MATLAB. In MALTAB, vectors are simply 1-row or 1-column matrices. In JULIA, vectors and matrices are two separate data types.

# Creating a Vector
# To create a vector, you can use square brackets [] and separate the elements by commas (NOT SPACE).
# Here's an example of creating a vector:
v = [1, 2, 3, 4, 5]

# Accessing Elements of a Vector
# You can access individual elements of a vector using indexing.
# In Julia, indexing starts from 1.
# Here's an example of accessing elements of a vector:
println(v[1])  # Output: 1
println(v[3])  # Output: 3

# Modifying Elements of a Vector
# You can modify individual elements of a vector by assigning a new value to the desired index.
# Here's an example of modifying elements of a vector:
v[2] = 10
println(v)  # Output: [1, 10, 3, 4, 5]

# Creating a Matrix
# To create a matrix, you can use square brackets [] and separate the rows by semicolons (;).
# Here's an example of creating a matrix:
m = [1 2 3; 4 5 6; 7 8 9]

# Accessing Elements of a Matrix
# You can access individual elements of a matrix using indexing.
# Similar to MATLAB, matrices are stored in column-major, with starting index 1.
# Here's an example of accessing elements of a matrix:
println(m[1, 2])  # Output: 2
println(m[3, 1])  # Output: 7
@show m[2:3,1]
@show m[2,:]

# Modifying Elements of a Matrix
# You can modify individual elements of a matrix by assigning a new value to the desired indices.
# Here's an example of modifying elements of a matrix:
m[2, 3] = 10
println(m)  # Output: [1 2 3; 4 5 10; 7 8 9]

# Creating view
# You can also create a view of a matrix, which is a reference to a subset of the original matrix.
msub = view(m, 2:3, 1:2)
@show msub[1, 1] 
# Because view is a direct reference to the memory address to where the original matrix is stored, modifying the view will also modify the original matrix.
msub[1, 1] = 100
@show m

# IMPORTANT: vectors are initialised by commas or semicolons separation. Row vectors are initialised by SPACE.
v = [1, 2, 3, 4, 5] # Column vector
rowMatrices = [1 2 3 4 5] # 1-Row matrix
# but you can always turn a 1-row matrix back into a vector
@show rowMatrices = vec(rowMatrices)

# Matrix Operations
# Matrix multiplication in Julia is done using the * operator, same as MATLAB.
A = [1 2; 3 4]
b = [1, 2]
c = A * b

# Matrix inversion is done using the \ operator, same as MATLAB.
x = A \ b

# One can also do element-wise operations using the . operator.
@show b.*[3,4]

# IMPORTANT: JULIA is stricter than MATLAB in specifying element-wise operations with functions (Broadcasting).
@show sin(b) # This will throw an error
@show sin.(b) # This will work

# IMPORTANT: JULIA is also stricter than MATLAB in initialising vector/matrices.
xvec = Vector{Float64}(undef, 5) # This will create a 5-element vector of zeros
for i in 1:5
    xvec[i] = i
end
@show xvec

## Tuples

# In Julia, a tuple is an ordered collection of elements, similar to an array.
# However, unlike arrays, tuples are immutable, which means their elements cannot be modified once defined.

# Creating a Tuple
# To create a tuple, you can use parentheses () and separate the elements by commas.
# Here's an example of creating a tuple:
t = (1,) # A tuple with one element
t = (1, 2, 3, 4, 5) # A tuple with multiple elements

# Accessing Elements of a Tuple
# You can access individual elements of a tuple using indexing.
# In Julia, indexing starts from 1.
# Here's an example of accessing elements of a tuple:
println(t[1])  # Output: 1
println(t[3])  # Output: 3

# Modifying Elements of a Tuple
# Since tuples are immutable, you cannot modify individual elements of a tuple.
# If you need to modify elements, you would need to create a new tuple with the desired changes.

# Tuple Unpacking
# In Julia, you can assign the elements of a tuple to separate variables using tuple unpacking.
# Here's an example of tuple unpacking:
a, b, c = (1, 2, 3)
println(a)  # Output: 1
println(b)  # Output: 2
println(c)  # Output: 3

println(t)          # Output: (1, 2, 3)

## Named Tuples
# Named tuples are similar to regular tuples, but with named fields.
# This allows you to access the elements of a named tuple using field names.

# Creating a Named Tuple
nt = (a=1,) # A named tuple with one element
nt = (a=1, b=2, c=3)

# Accessing Elements of a Named Tuple
# You can access individual elements of a named tuple using field names.
# Here's an example of accessing elements of a named tuple:
println(nt.a)  # Output: 1

## IMPORTANT: tuples are immutable, which means their elements cannot be modified once defined.
nt.a = 10 # This will throw an error
# but you can always redefine a new tuple with the same variable named
nt = (a=10, b=2, c=3, d=4)

# Symbols in tuples
# The name of the field in a named tuple is a "symbol", which one can think of as immutable names to certain internal variable.
# Symbols are created by prefixing a variable with a colon (:). For example, you can call field "a" in the named tuple nt as :a
println(nt[:a]) # Output: 10 

## Functions

# Julia allows you to define your own functions using the function keyword.

# Here's an example of a function that calculates the square of a number:
function square(x)
    return x^2
end

# You can call the function like this:
@show square(5)

# You can also define functions with named optional arguments. Arguments after a semicolon are optional, and you can specify default values for them.
function power(x; n=2)
    return x^n
end

@show power(5)
@show power(5, n=3)

# Mutating Functions
# Usually, we recommend poeple to write functions in Julia that does not modify their arguments.
# However, if you want to modify the argument, you can use the keyword "!" at the end of the function name to notify the user that the function is mutating.

function mutate!(x,y)
    x[1] = 10 
    return x.+y
end

aa=[2, 3]
bb=[3, 5]

@show mutate!(aa,bb)
@show aa

# IMPORTANT: However, note that an assignment to a variable always define a new variable, and will therefore not mutate the original variable.

function NotMutate!(x,y)
    x = 10 
    return x.+y
end
aaa = 2
bbb = 3
@show NotMutate!(aaa,bbb)
@show aaa

## Macros
# Macros is a feature in JULIA that is similar to Macros in C. They are lower-level functions that we don't usually recommend for beginners to use.
# But a few macros are quite handy

# @show
# The @show macro is used to print the expression and its value. Since by default, code in script is not printed, this is a good way to show evaluted values in the middle of the code.
@show x

# @time
# The @time macro is used to measure the time taken to execute a line/block of code.
@time 1+2

## Ranges and Plotting
# In Julia, you can create ranges using the colon operator (:).
# A range is a sequence of numbers that are evenly spaced within a specified interval.
# The syntax is the same as MATLAB.
# They are practically the same as arrays, but they are more efficient.
xmesh = 0:0.1:π

# To plot things, you can use the Plots package.
using Plots
plotlyjs() # You can choose a backend. plotlyjs is a nice feature-rich backend.
plot(xmesh, sin.(xmesh); label="sin")
plot!(xmesh, cos.(xmesh); label="cos") # To overlap more plots, use plot! instead of plot