# PRaCTES session 7: high performance computing

"High performance computing" usually means writing and running parallel programs on large numbers of processors. I'm interpreting the scope of this session a little more broadly, though. The first half will focus on understanding the performance of serial programs, with an emphasis on why some programming languages produce code that is dramatically (i.e., several orders of magnitude) faster tha others. The second half will focus on understanding how parallel programs decrease program execution times and provide some practing running and interpreting benchmarks to measure the speedup produced by parallelization.

## Part 1: programming languages, compilers, and serial performance

We'll explore the performance of serial programs with a toy program based on the Collatz conjecture, which proposes that the iterative sequence

<a href="https://www.codecogs.com/eqnedit.php?latex=n&space;\rightarrow&space;n/2&space;\;\;\textrm{&space;(even&space;}&space;n&space;\textrm{)}&space;\\&space;n&space;\rightarrow&space;3n&plus;1&space;\;\;\textrm{&space;(odd&space;}&space;n&space;\textrm{)}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?n&space;\rightarrow&space;n/2&space;\;\;\textrm{&space;(even&space;}&space;n&space;\textrm{)}&space;\\&space;n&space;\rightarrow&space;3n&plus;1&space;\;\;\textrm{&space;(odd&space;}&space;n&space;\textrm{)}" title="n \rightarrow n/2 \;\;\textrm{ (even } n \textrm{)} \\ n \rightarrow 3n+1 \;\;\textrm{ (odd } n \textrm{)}" /></a>

always reaches 1. For example, the starting value 6 produces a sequence (6, 3, 10, 5, 16, 8, 4, 2, 1) that eventually reaches 1 after 8 iterations. Our toy problem is the following: given some number *N*, what starting value between 1 and *N* takes the most iterations to reach 1?

The ``collatz_comparison`` folder contains a solutions to this problem written in four programming languagues: Python, C, Julia, and Java. Although all four programs produce the same answers, they run at very different speeds, and we'll try to understand why---but first, let's talk about how computers execute code.

### How computers execute code

Computers contain a large bank of memory ("random access memory", or RAM), which can typically holds several GB of data, and one or more processors. The processors themselves contain several much smaller banks of memory (called "registers"), which hold just 64 bits of data on most modern computers, and circuits of transistors that can perform basic arithmetic and logic operations on values stored in the registers. Ultimately, all programs run by executing sequences of instructions ("machine code") that manipulate data in registers and move data between registers and RAM.

For example, a function that adds two integers
```C
int add(int a, int b) {
   return a + b;
}
```
might be converted to a sequence of instructions that looks something like
```
load a, r1          (load a from memory into register 1)
load b, r2          (load b from memory into register 2)
add r1, r2, r3      (add the data in registers 1 and 2 and store the result in register 3)
store r3, retval    (store the data in register 3 in a special "return value" location in RAM)
```

Different types (or "architectures") of processors use different instructions sets, so the same program will be converted to different machine code if it runs on e.g. an Intel processor (which uses the x86 instruction set) and an IBM processor (which uses the Power instruction set). Generally speaking, programs will run more quickly if they require a relatively short sequence of instructions.

### Language 1: Python

Python programs are never converted directly to machine code. Instead, they are translated into machine-independent code objects that run on the "Python virtual machine"---which is itself a program, often written in C. Because the translation process doesn't depend on the instruction set used by the computer hardware, the Python virtual machine makes it much easier to run Python on a variety of processor architectures. However, the extra layer of software between the program and the hardware introduces some overhead that slows program execution.

Start a Python REPL and measure the time to solution for *N* = 1,000,000 by running
```python
>>> from collatz import *
>>> import timeit
>>> def benchmark():
...    find_longest_sequence(1000000)
>>> timeit.timeit(benchmark, number = 1)
```
This takes 1-2 minutes on my laptop and tells me that the longest sequence starts from 837,799 and reaches 1 after 524 iterations.

### Language 2: C

C is a compiled language. This means that C programs are first converted all the way to machine code by running them through a program called a compiler, and the machine code is then executed directly by the processor. This avoids virtual machine overhead and produces programs that execute quickly. However, it also requires developing and maintaining many versions of the C compiler, since a different version is required for every processor architecture.

From the command line, compile the C program and measure the time to solution for *N* = 1,000,000 by running
```bash
$ gcc collatz.c -o collatz
$ time ./collatz 1000000
```
This takes around 1 second on my laptop---much faster than the Python program!

These two experiments provide the basis for a hypothesis:

### Programs run faster if they compile to machine code

We are scientists, so let's test this hypothesis by collecting more data. Specifically, let's measure the time to solution in two other programming languages: Java (which runs on a virtual machine) and Julia (which is compiled to machine code).

### Language 3: Java

From the command line, compile the Java program (to machine-independent code objects, not machine code!) and run the program (with ``-Xint`` to disable just-in-time compilation to machine code).
```bash
$ javac Collatz.java
$ time java -Xint Collatz 1000000
```
On my laptop, this takes about 20 seconds.

### Language 4: Julia

Start a Julia REPL and measure the time to solution by running
```julia
julia> include("collatz.jl")
julia> find_longest_sequence(1)
julia> @time find_longest_sequence(1000000)
```
On my laptop, this also takes about 20 seconds. (The first call to ``find_longest_sequence`` is to trigger Julia's just-in-time compiler so that our measurement doesn't include the time needed for the compiler to run.)

This new data suggests that our hypothesis is incomplete. Like Python, Java runs on a virtual machine, but it runs much more quickly. Additionally, Julia runs much slower than C (and comparable to Java) even though it compiles to machine code.

## Exercise 1: 

## Part 2: parallel computing
