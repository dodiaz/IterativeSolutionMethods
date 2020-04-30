# Implementation of different iterative methods to solve poisson's equation

Each of the 5 methods are implemented in one script called "IterativeSolvers":
1. Point Jacobi
2. Gauss-Seidel
3. Successive over-relaxation
4. Conjugate gradient
5. Multigrid

There is a variable called "method" that needs to be changed when running the code to test any of the iterative methods. The possible methods are: "PJ", "GS", "SOR", "CG", or "MG".

I've altered the function "print_current_data" so that it prints laplace_phi, f, phi, and error but I haven't tested it out so I'm not sure if it works yet.