## ----------------------------------------------------------------- 
# load packages
# ------------------------------------------------------------------
using CairoMakie
using Zygote
using Optimization
using OptimizationOptimJL
using LinearAlgebra
using Base.Threads
using FileIO, JLD2
using LineSearches
using ProgressMeter
using LazySets
using LazySets: Interval
using ArgParse
using Dates
using Base.Iterators

## --------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
# Base period (in days) M of the types. This is multiplied by a factor L to 
# get the actual period.
baseM = 30                  

# Multiplying factor for the period L ranges from 1 to maxL. With maxL = 2, 
# there will be types with period baseM and 2*baseM. 
maxL = 1                    

# Frequency of the environmental payoff factor within a cycle with length baseM
freq = 2                    

# Range of values of the selection strength w
w_begin = 0.75
w_end = 0.75
w_num = 1

# Range of values of the shape parameter μ
μ_begin = 0.01
μ_end = 1.0
μ_num = 100

# List of values of the spread parameter δ
# A vector (list) of integers is applicable when maxL = 1. 
# When maxL > 1, each element of δ_list should be a vector (list) of maxL 
# elements. Each corresponds to the spread parameter for the appropriate 
# type of baseM*L period.
δ_list = collect(1:15)

# The phase of the tides, which factors in the environmental fitness factor
θ = 8.25

# The parameter controlling the non-linearity of the density dependence 
# of reproductive benefits
h = 0.

# The solver uses a local algorithm (Newton's method) to find the solution,
# as it works the best for such problems and in finding all solutions. 
# This algorithm needs initial guesses in order to converge to the solution. 
# numPoints_1 is the default number of points per side  of the M-simplex to 
# use as the default number of guess points, which are all uniformly spaced. 
numPoints_1 = 1

# There are some regions of parameter space which the solver finds difficult
# to solve (determined visually through discontinuities in roughness in the 
# solution plots). For regions, a higher number of initial guesses is required
# to find the true number of solutions. switch_guess_points indicates whether
# there are problem areas in the parameter space. problem_area designates the
# approximate region of the parameter space which needs to be treated more 
# carefully. Some options are provided which work well for the parameters 
# of the figures. numPoints_2 is the number of guess points per side of the 
# simplex for parameters in the problem area. numPoints_L designates an even 
# larger number of guess points for especially problematic conditions.
switch_guess_points = true
problem_area_options = (
    [Interval(0.75, 1.0), Interval(0.0, 2.0), Interval(0.0, 8.0)],
    [Interval(0.75, 1.0), Interval(0.0, 2.0), Interval(0.0, 8.0) × Interval(6.0, 15.0)],
    [Interval(0.6, 1.0), Interval(0.0, 2.0), Interval(0.0, 20.0)],
)  #μ, w, δs (cartesian product of intervals for individual δs)
problem_area_number = 1

numPoints_2 = 2
numPoints_L = 3

# As the solving process proceeds to find solutions to more parameter values,
# The solutions of neighbouring parameters are used as initial guesses for 
# the local solving method. This makes the solving faster and more robust.
# neighbours sets the number of closest neighbours to use solutions from.
neighbours = 10

# Tolerance of the solver. When iterative changes to solution are smaller
# than this value, the solving method stops and returns the solution.
solTolerance = 1e-12

# Tolerance to consider two solutions obtained from different initial
# guesses are unique.
uniquenessTolerance = 1e-6

# Name of the file and the path where the results are to be saved.
savefile = "M30_solutions"
pathname = "Data/"

## ----------------------------------------------------------------- 
# function definitions
# ------------------------------------------------------------------
# Non-linear density dependent reproductive factor of payoff
PayoffFunction(x, h) = x*(1+h)/(1+h*x)
# First derivative of the reproductive payoff factor
PayoffFunctionPrime(x, h) = (1+h)/(1+h*x) - (h*x)*(1+h)/(1+h*x)^2

# Augmented equations to be solved for. Plugging in the type frequencies,
# payoff vector, mutation matrix, it returns the error consisting of:
# 1. Squared sum of the RHS of the equation
# 2. Normalization error to keep the frequency normalized
# 3. Bounding error to keep the frequency within valid bounds (0, 1)
# When passing in a valid solution, this function returns 0.
function equations(X, (payoffVector, mutationMatrix, h) )
    ϕ = sum(payoffVector .* X .* PayoffFunction.(X, h)) 
    eqs = transpose(mutationMatrix)*(payoffVector .* X .* PayoffFunction.(X, h)) - ϕ*X
    return sum(eqs .^ 2) + ( sum(X)-1 )^2  + sum( (X .< 0) .* X.^2 ) + sum( (X .> 1) .* (X .- 1).^2 )
end

# Jacobian matrix of the above function for the Newton's method solver 
# and to determine stability. 
function jacobian(X, payoffVector, mutationMatrix, h, M)
    ϕ = sum(payoffVector .* X .* PayoffFunction.(X, h)) 
    a = payoffVector
    Q = mutationMatrix
    J = zeros(M, M)
    for i in 1:M
        for j in 1:M
            J[i,j] = a[j]*(PayoffFunction(X[j], h) + X[j]*PayoffFunctionPrime(X[j], h))*(Q[j,i] - X[i]) - a[M]*(PayoffFunction(X[M], h) + X[M]*PayoffFunctionPrime(X[M], h))*(Q[M,i] - X[i])
            if i==j
                J[i,i] -= ϕ
            end
        end
    end
    return J
end

# Determines if a solution is stable using the above Jacobian function.
# A solution is stable if all the eigenvalues of the Jacobian matrix
# have negative real part. This function returns the maximum real 
# part of the eigenvalues.
function stability(X, payoffVector, mutationMatrix, h, M)
    J = jacobian(X, payoffVector, mutationMatrix, h, M)
    if all(isfinite, J)
        return maximum(real.(eigvals(J)))
    else
        return 1000
    end
end

# Recursive function to find a grid of equally spaced
# guesspoints within a simplex of dimension `dim`.
function SimplexPointsLower(dim, pointsPerSide)
    if dim == 1
        return collect(0:1/pointsPerSide:1)
    end
    S_lower = SimplexPointsLower(dim-1, pointsPerSide)
    S = Float64[]
    for i in axes(S_lower, 1)
        point = S_lower[i,:]
        for j in 0:1/pointsPerSide:(1-sum(point))
            if size(S,1)==0
                S = [S; [point; j]]
            else
                S = [S [point; j]]
            end
        end
    end
    return S'
end

# The above function finds a relatively minimal representation of the
# full set of guesspoints. That minimal representations is used to 
# calculate the remaining set of guesspoints.
function SimplexPointsFull(dim, pointsPerSide)
    lowerSimplex = SimplexPointsLower(dim-1, pointsPerSide)
    fullSimplex = zeros(Float64, size(lowerSimplex, 1), dim)
    fullSimplex[:, 1:end-1] .= lowerSimplex
    for i in 1:size(lowerSimplex, 1)
        fullSimplex[i, end] = 1 - sum(fullSimplex[i, 1:end-1])
    end
    return fullSimplex
end

# This function calculates the grid of guesspoints when there are types
# with different periodicities involved in the equations.
function MultipleSimplexPoints(baseM, maxL, pointsPerSide, mixing=false)
    fullM = baseM*maxL*(maxL+1)÷2
    if mixing == true
        return SimplexPointsFull(fullM, pointsPerSide)
    end
    all_points_compact = [SimplexPointsFull(baseM*period, pointsPerSide) for period in 1:maxL]
    all_points = [zeros(Float64, size(m,1), fullM) for m in all_points_compact]
    index = 0
    for L in 1:maxL
        all_points[L][:, 1+index:index + baseM*L] = all_points_compact[L]
        index += baseM*L
    end
    all_points = reduce(vcat, all_points)
    return all_points
end

# This function pushes all the guesspoints that are on any boundary of the
# simplex slightly into the interior so that they are valid guesses for 
# Newton's method.
function interiorize(guessPoints)
    S = copy(guessPoints)
    l,M = size(S)
    if M == 1
        S[1] += 1e-6
        S[end] -= 1e-6
    else
        for i in 1:l
            for j in 1:M
                if S[i,j] == 0.
                    S[i,j] += 1e-6
                elseif S[i,j] == 1.
                    S[i,j] -= 1e-6
                end
            end
        end
    end
    return S
end

# This function removes duplicate solutions from a list of solutions
# which differ by more than the tolerance `tol` in their norm-2 difference.
function unduplicate(sols, tol=1e-6)
    indicesToUse = Int64[]
    indicesToReject = Int64[]
    for i in 1:size(sols,2)
        if !(i in indicesToReject)
            indicesToUse = [indicesToUse; i]
            for j in i+1:size(sols,2)
                if !(j in indicesToReject) && sum((sols[:,i] .- sols[:,j]).^2) < tol
                    indicesToReject = [indicesToReject; j]
                end
            end
        end
    end
    sols[:,indicesToUse]
end

# This function solves the equations when given all parameters and guess points
function solveEqs(guessPoints, payoffVector, order, mutationMatrix, h, M, soltol=1e-6, uniqtol=1e-6)
    num = size(guessPoints, 2)
    sols = zeros(M,num)
    optf = OptimizationFunction(equations, AutoZygote())

    # multithreaded solving with load distribution
    chunks = Iterators.partition(axes(guessPoints, 2),max(num÷Threads.nthreads(), 1))
    tasks = map(chunks) do cols
        @spawn begin
            for col in cols
                prob = OptimizationProblem(optf, guessPoints[:,col], (payoffVector, mutationMatrix, h))
                try
                    sols[:,col] = solve(prob, Optim.Newton(; alphaguess = InitialHagerZhang()))
                catch err
                    println(err)
                    sols[1,col] = -10
                end 
            end
        end
    end
    wait.(tasks)

    # removing the solutions which did not converge
    sols = sols[:, sols[1,:] .!= -10]
    # removing the unstable solutions
    st = [stability(sols[:,i], payoffVector, mutationMatrix, h, M) for i in axes(sols, 2)]
    sols = sols[:, findall(st .< 0)]
    # removing repeating solutions
    sols = unduplicate(sols, uniqtol)
    # removing solutions that are not good enough and keeping only the best M
    errors = [equations(sols[:,i], (payoffVector, mutationMatrix, h)) for i in axes(sols, 2)]
    pass_indices = findall(errors .< soltol)
    if size(pass_indices, 1) > M
        sols = sols[:, sortperm(errors)[1:M]]
    else
        sols = sols[:, pass_indices]
    end
    # rearranging the solutions in order of decreasing fitness
    if size(sols, 2) > 1
        numsols = size(sols, 2)
        maxIndices = [findmax(sols[:,i])[2] for i in 1:numsols]
        rearr_indices = [findfirst(order[i] .== maxIndices) for i in axes(order, 1)]
        rearr_indices = rearr_indices[rearr_indices .!= nothing]
        sols = sols[:, rearr_indices]
    end
    return sols
end

# Environmental payoff factor vector with selection strength w, tidal phase θ,
# cycle length M and frequency within that `freq`
function cosinePayoffVector(w, θ, M, freq = 1) 
    A = 1 .+ w * cos.(2*freq*pi/M*(collect(1:M) .- θ))
    A[A .< 0.] .= 0.
    return A
end

# Environmental payoff vector when types with multiple periodicities are involved.
function multipleCosinePayoffVector(w, θ, baseM, maxL, freq = 1)
    payoffVector = reduce(vcat, [cosinePayoffVector(w, θ, baseM*period, freq*period) for period in 1:maxL])
    return payoffVector
end

# This function is used to order obtained solutions so that they are in descending 
# order of fitness and subsequent ones differ incrementally in phase.
function multipleOrder(θ, baseM, maxPeriod)
    order = zeros(Int64, baseM*maxPeriod*(maxPeriod+1)÷2)
    for i in 1:maxPeriod
        start = i*(i-1)÷2*baseM + 1
        fin = i*(i+1)÷2*baseM
        maxPayoff = cosinePayoffVector(1, θ, baseM*i, i)
        order[start:fin] = sortperm(maxPayoff)[end:-1:1] .+ i*(i-1)÷2*baseM
    end
    return order
end

# This function creates a banded mutation matrix with shape parameter μ,
# spread parameter δ and cycle length M
function bandedMutationMatrix(μ, δ, M)
    if δ>M/2-1
        mutationVector = [1-μ+μ/M; ones(M-1)*μ/M]
    else
        mutationVector = [1-μ+μ/(2δ+1); ones(δ)*μ/(2δ+1); zeros(M-2δ-1); ones(δ)*μ/(2δ+1)]
    end
    mutationMatrix = [mutationVector[abs(j-i)%M+1] for i in 1:M, j in 1:M]
    return mutationMatrix
end

# This function creates the mutation matrix when types of different periodicities 
# are involved.
function multipleBandedMutationMatrix(μ, δs, baseM, maxL)
    fullM = baseM*maxL*(maxL+1)÷2
    mutationMatrix = zeros(Float64, fullM, fullM)
    index = 0
    for period in 1:maxL
        mutationMatrix[1+index:index + baseM*period, 1+index:index + baseM*period] .= bandedMutationMatrix(μ,δs[period],baseM*period)
        index += baseM*period
    end
    return mutationMatrix
end

# This function finds the closes N neighbours to a point in parameter space
# given a grid of points. This is required to use the neighbouring solutions as guesses.
function closestNNeighbours((x,y,z), presenceMatrix, number)
    numx,numy,numz = size(presenceMatrix)
    distance_list = [[(i-x)^2+(j-y)^2+(k-z)^2;i;j;k] for i in 1:numx, j in 1:numy, k in 1:numz if presenceMatrix[i,j,k] && ~(i==x && j==y && k==z)]
    if size(distance_list,1) == 0
        return [;;;]
    end
    distance_list = permutedims(reduce(hcat,distance_list),(2,1))
    order = sortperm(distance_list[:,1])
    distance_list = distance_list[order,:]
    closest = [CartesianIndex(distance_list[i,2],distance_list[i,3],distance_list[i,4]) for i in 1:min(number,size(distance_list,1))]
    return closest
end

# When there are types with different periodicities in the same system,
# There are types which emerge on the same day. This function condenes 
# the full solution into one where each day of the cycle is represented 
# only ones. All types, regardless of what their periodicity is are 
# added up on the day which they emerge. 
function compactSolutions(sols, baseM, maxL, uniq_tol)
    comL = lcm(1:maxL)
    comsols = zeros(comL*baseM, size(sols,2))
    for i in axes(sols, 2)
        sol = sols[:,i]
        for period in 1:maxL
            st = 1+period*(period-1)÷2*baseM
            fin = st - 1 + baseM*period
            comsols[:,i] .+= repeat(sol[st:fin], comL÷period)./(comL÷period)
        end
    end
    comsols = unduplicate(comsols, uniq_tol)
    return comsols
end

# This function wrapps the solver function and takes care of the problem cases
# and tries to solve them with higher guesspoints.
function solverWrapper!(solutions, comsolutions, solved, pars, is_problem, guessPoints_1, guessPoints_2, guessPoints_L, neighbours, order, h, fullM, solTolerance, uniquenessTolerance, pbar, (i,j,k))
    if !is_problem
        if ~any(solved)
            new_guessPoints = guessPoints_L
        else
            new_guessPoints = hcat(guessPoints_1, solutions[closestNNeighbours((i,j,k), solved,neighbours)]...)
            new_guessPoints = unduplicate(new_guessPoints, uniquenessTolerance)
        end
        solutions[i,j,k] = solveEqs(new_guessPoints, pars[1], order, pars[2], h, fullM, solTolerance, uniquenessTolerance)
    else
        if ~any(solved)
            new_guessPoints = guessPoints_L
        else
            new_guessPoints = hcat(guessPoints_2, solutions[closestNNeighbours((i,j,k), solved, neighbours)]...)
            new_guessPoints = unduplicate(new_guessPoints, uniquenessTolerance)
        end
        solutions[i,j,k] = solveEqs(new_guessPoints, pars[1], order, pars[2], h, fullM, solTolerance, uniquenessTolerance)
    end
    comsolutions[i,j,k] = compactSolutions(solutions[i,j,k], baseM, maxL, uniquenessTolerance)
    solved[i,j,k] = true
    next!(pbar)
    GC.safepoint()
    return nothing
end

# This function finds all solutions of each points of a grid of points in 
# the parameter space
function gridSolve(pars, guessPoints_1, guessPoints_2, guessPoints_L, is_problem, neighbours, order, w_num, μ_num, δ_num, h, baseM, maxL, solTolerance, uniquenessTolerance, pathname, savefile)
    fullM = baseM*maxL*(maxL+1)÷2
    solutions = Array{Matrix}(undef, w_num, μ_num, δ_num)
    comsolutions = Array{Matrix}(undef, w_num, μ_num, δ_num)
    solved = zeros(Bool, w_num, μ_num, δ_num)
    tasks = Array{Task}(undef, w_num, μ_num, δ_num)
    open(pathname*savefile*"_prog.txt","a") do io
        pbar = Progress(w_num*μ_num*δ_num; output = io)
        update!(pbar,0)
        for i in 1:w_num
            for j in 1:μ_num
                for k in 1:δ_num
                    if ~is_problem[i,j,k]
                        tasks[i,j,k] = @spawn solverWrapper!(solutions, comsolutions, solved, pars[i,j,k], is_problem[i,j,k], guessPoints_1, guessPoints_2, guessPoints_L, neighbours, order, h, fullM, solTolerance, uniquenessTolerance, pbar, (i,j,k))
                    else
                        tasks[i,j,k] = @spawn 1
                    end
                end
            end
        end
        wait.(tasks)
        for i in 1:w_num
            for j in 1:μ_num
                for k in 1:δ_num
                    if is_problem[i,j,k]
                        tasks[i,j,k] = @spawn solverWrapper!(solutions, comsolutions, solved, pars[i,j,k], is_problem[i,j,k], guessPoints_1, guessPoints_2, guessPoints_L, neighbours, order, h, fullM, solTolerance, uniquenessTolerance, pbar, (i,j,k))
                    end
                end
            end
        end
        wait.(tasks)
        finish!(pbar)
    end
    return solutions, comsolutions
end

# This function parses parameters provided as command line flags 
function parseCommandLine()
    s = ArgParseSettings(add_help = false)
    @add_arg_table(s,
        "--baseM", arg_type = Int, default = baseM,
        "--maxL", arg_type = Int, default = maxL,
        "--freq", arg_type = Int, default = freq,
        "--w_begin", arg_type = Float64, default = w_begin,
        "--w_end", arg_type = Float64, default = w_end,
        "--w_num", arg_type = Int, default = w_num,
        "--mu_begin", arg_type = Float64, dest_name = "μ_begin", default = μ_begin,
        "--mu_end", arg_type = Float64, dest_name = "μ_end", default = μ_end,
        "--mu_num", arg_type = Int, dest_name = "μ_num", default = μ_num,
        "--switch_guess_points", arg_type = Bool, default = switch_guess_points,
        "--theta", arg_type = Float64, dest_name = "θ", default = θ,
        "--h", arg_type = Float64, default = h,
        "--numPoints_1", arg_type = Int, default = numPoints_1,
        "--problem_area_number", arg_type = Int, default = problem_area_number,
        "--numPoints_2", arg_type = Int, default = numPoints_2,
        "--numPoints_L", arg_type = Int, default = numPoints_L,
        "--neighbours", arg_type = Int, default = neighbours,
        "--solTolerance", arg_type = Float64, default = solTolerance,
        "--uniquenessTolerance", arg_type = Float64, default = uniquenessTolerance,
        "--savefile", arg_type = String, default = savefile,
        "--pathname", arg_type = String, default = pathname,
        "--evalCode", arg_type = String, default = "",
    )
    return parse_args(s)
end

## ----------------------------------------------------------------- 
# Initial setup
# ------------------------------------------------------------------
parsed_args = parseCommandLine()
for (arg, val) in parsed_args
    @eval $(Symbol(arg)) = $val
end

eval(Meta.parse(evalCode))
println()

fullM = baseM*maxL*(maxL+1)÷2
comM = baseM*lcm(1:maxL)

δ_num = length(δ_list)
if all([μ_num, w_num, δ_num] .> 1)
    error("Only two out of μ, w and δ can be a range.")
end

if w_num == 1
    w_range = [w_begin]
else
    w_range = range(w_begin,w_end,w_num)
end

if μ_num == 1
    μ_range = [μ_begin]
else
    μ_range = range(μ_begin,μ_end,μ_num)
end

problem_area = problem_area_options[problem_area_number]

# Create the grid of guess points from the number per side
guessPoints_1 = MultipleSimplexPoints(baseM,maxL,numPoints_1);
guessPoints_1 = interiorize(guessPoints_1)';
guessPoints_2 = MultipleSimplexPoints(baseM,maxL,numPoints_2);
guessPoints_2 = interiorize(guessPoints_2)';
guessPoints_L = MultipleSimplexPoints(baseM,maxL,numPoints_L);
guessPoints_L = interiorize(guessPoints_L)';
pars = [(multipleCosinePayoffVector(w,θ,baseM,maxL,freq), multipleBandedMutationMatrix(μ,δs,baseM,maxL)) for w in w_range, μ in μ_range, δs in δ_list];

# Create the vector that is used to order solutions appropriately
ordering = multipleOrder(θ,baseM,maxL)
if switch_guess_points
    is_problem = [[μ] ∈ problem_area[1] && [w] ∈ problem_area[2] && δs ∈ problem_area[3] for w in w_range, μ in μ_range, δs in δ_list]
else
    is_problem = zeros(Bool,w_num,μ_num,δ_num)
end

try
    Base.Filesystem.mkdir(pathname)
catch
end

# Print out the parameters selected. Progress of the calcualations is printed to the same file
open(pathname*savefile*"_prog.txt","w") do file
    write(file,string(now()),"\n")
    for key in ["baseM","freq","maxL","w_begin","w_end","w_num","μ_begin","μ_end","μ_num","θ","δ_list","h","switch_guess_points","problem_area","numPoints_1","numPoints_2","numPoints_L","numPoints_L","neighbours","solTolerance","uniquenessTolerance","savefile","pathname"]
        value = string(eval(Symbol(key)))
        write(file,key," => ",value,"\n")
    end
end

## ----------------------------------------------------------------- 
# Solving the equations
# ------------------------------------------------------------------
time1 = time()
fullsolutions, comsolutions = gridSolve(pars, guessPoints_1,guessPoints_2, guessPoints_L,is_problem, neighbours,ordering,w_num,μ_num,δ_num,h,baseM,maxL,solTolerance,uniquenessTolerance,pathname,savefile);

time2 = time()

numsols = zeros(Int64,w_num,μ_num,δ_num)

# repackaging of the solutions
solutions = zeros(Float64,comM,fullM,w_num,μ_num,δ_num) # solution, solution idx, w, μ, δ
for i in 1:w_num
    for j in 1:μ_num
        for k in 1:δ_num
            ns = size(comsolutions[i,j,k],2)
            numsols[i,j,k] = ns
            for u in 1:ns
                solutions[:,u,i,j,k] .= comsolutions[i,j,k][:,u]
            end
        end
    end
end

println("Elapsed Time: ", round(time2-time1,digits=2)," s") 

## ----------------------------------------------------------------- 
# Saving the solutions and parameters
# ------------------------------------------------------------------

jldopen(pathname*savefile*".jld2","w") do file
    for key in ["baseM","maxL","freq","w_begin","w_end","w_num","μ_begin","μ_end","μ_num","θ","δ_list","δ_num","h","numPoints_1","numPoints_2","numPoints_L","neighbours","switch_guess_points","problem_area","solTolerance","uniquenessTolerance","savefile","pathname","solutions","numsols","fullsolutions","fullM"]
        value = eval(Symbol(key))
        write(file,key,value)
    end
end