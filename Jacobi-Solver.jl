using DataFrames, Statistics, DelimitedFiles, CSV, LinearAlgebra, Random, Distributions
using Pkg
Pkg.update()
#pseudo data generation 5x5 matrix
b = randn(5)
X0 = zeros(5)
A = [1 2 0 1 3; 0 1 2 1 0; 1 2 1 2 2; 1 4 2 2 1; 2 1 2 0 1]

#jacobi solver
# """function for jacobi solver
# b: is the response vector
# A: is the design matrix for intercept + slope
# iter: is the number of the iterations in the jacobi
# Tol: Tolerance
# X0: Initial guess for k variables
# """
function jacobi(A,b,N,x,TOL)
    """Solves the equation Ax=b via the Jacobi iterative method."""
    # Create an initial guess if needed
    # if x
    #     x = zeros(size(A, 2))
    # end
    # Create a vector of the diagonal elements of A
    # and subtract them from A
    D = diag(A)
    R = A .- D
    println(size(R))

    # Iterate for N times
    for i in 1:N
        X0 = x
        x = inv(Diagonal(D))*(b .- R*x)
        if abs.((x .- X0)'*(x .- X0)) .< TOL
            return round.(x, digits=2)
        end
    end
    return round.(x, digits=2)
end


#small genotype matrix 14000 x 100
geno = readdlm("../SSBR-JWAS-Implementation\\SSBR-vs-BR-Example\\data\\Input-Data\\Geno_Small.txt")
println(size(geno))
ped = readdlm("../SSBR-JWAS-Implementation\\SSBR-vs-BR-Example\\data\\Input-Data\\Pedigree_XSim_Raw1.txt")

#get A matrix
s=ped[:,2]
println(s)
d=ped[:,3]
n = length(s)
s=(s .== 0)*n .+s
d=(d .== 0)*n .+d;
A = zeros(n,n);
println(s[2])
for i in 1:n
    A[i,i] = 1 + A[Int64(s[i]), Int64(d[i])]/2
    for j in (i+1):n
        A[i,j] = ( A[i, Int64(s[j])] + A[i, Int64(d[j])] ) / 2
        A[j,i] = A[i,j]
    end
end
println(A[1, 1])


#Define genetic and residual variances
vara = Diagonal(repeat([2.0], size(geno, 2)))
vare = Diagonal(repeat([1.0], size(geno, 1)))

# Define random phenotype of animals
y = randn(size(geno, 1))

#Construct MME
lhs = (geno'*inv(vare)*geno + inv(vara))
rhs = geno'*inv(vare)*y

#Solve MME using jacobi solver
X0 = randn(size(geno, 2))
alpha = jacobi(lhs, rhs, 50, X0, 0.00002)
println(alpha[1])
#Get y pred and compare to y observed
y_pred = geno*alpha .+ diag(vare)

#Condition of MME matrix
cond_mme = cond(lhs)
#calculate correlation between y_pred and y_obs
accuracy = cor(y, y_pred)
