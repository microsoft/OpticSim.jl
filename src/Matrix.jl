# Keep this code as a library for future use. The JacobiSVD is significantly faster than the SVD that comes standard with Julia.

function ATA(A::AbstractArray{T,2}, j::Int, k::Int) where {T<:Real} # need AbstractArray instead of Array to handle Transpose, which is an AbstractArray but not an Array
    sum = zero(T)
    nrows, _ = size(A)
    @simd for i in 1:nrows
        sum += A[i, j] * A[i, k]
    end
    return sum
end

@inline epsilon(T::Real) = 4 * eps(T)

"this version will always perform iterations to high relative accuracy"
GivensRotation(ajj::T, ajk::T, akk::T) where {T<:Real} = GivensRotation(false, ajj, ajk, akk)

function GivensRotation(lowRelativeAccuracy::Bool, ajj::T, ajk::T, akk::T) where {T<:Real}
    local sinTheta, cosTheta, needToRotate

    product = ajj * akk # default value that may be overriden if computing small singular values to low relative accuracy

    if lowRelativeAccuracy
        println("using low relative accuracy")
        # less stringent accuracy requirement than requiring that the rotation matrices not change the diagonal elements of the matrix at all. The result is that small singular values may be computed with poor 
        # relative accuracy, possibly having few significant digits of accuracy. For most applications computing the small singular values to high accuracy is not important so this is probably the way you would 
        # want to compute most of the time.
        max = max(abs(ajj), abs(akk))
        min = min(abs(ajj), abs(akk))
        if min == 0 || max / min > 2^52  # this will cause the test "if (abs(ajk) >= eps * sqrt(product))" to fail with fewer iterations which makes the code significantly more efficient. This needs to be rewritten in terms of machine epsilon. This number pow(2,52) is correct only for Float64. 
            product = max
        end
    end

    if ajk != zero(T) && abs(ajk) >= eps(T) * sqrt(product)
        local c, s, t, tau

        tau = (ajj - akk) / (2 * ajk)
        t = sign(tau) / (abs(tau) + sqrt(1 + tau * tau))
        c = 1 / sqrt(1 + t * t)
        s = c * t
        sinTheta = s
        cosTheta = c

        # NOTE: commented out this code because it greatly increased the residual Ax-b. Not sure how to include this without messing up accuracy in general. Leave it out for now.

        # when the matrix is singular it can be the case that the diagonal elements ajj, akk go to zero faster than ajk. This means the test "if (abs(ajk) >= eps * sqrt(ajj*akk))"
        # returns true until the diagonal and off diagonal elements become so small that the rotation matrix has no effect on any of the elements in the matrix (the rotation matrix becomes almost exactly
        # equal to the identity). If this happens the JacobiSVD function 
        # will never stop iterating and you will get a max iterations exceeded exception. The following test detects when the diagonal elements are no longer being changed by rotation
        # and causes the iteration to stop if this is the case.
        # if abs((c * ajj - s * ajk) - ajj) < eps(T) && abs((s * ajk + c * akk) - akk) < eps(T)
        #     # c and s values are such that none of the diagonal elements in the matrix will change if the rotation is performed so don't rotate. More iterations than necessary since have to wait till diagonal and
        #     # off diagonal elements are reduced to values so small they are beyond the range of whatever precision float you are using. This is probably necessary if want high relative accuracy in small singular values but
        #     # there are many cases when this is not a requirement.
        #     # Have to use eps(T) instead of a test for equality because C# is apparently not using 64 bit IEEE double precision arithmetic. Instead it appears to be using the Intel extended precision 
        #     # floating point. It takes many more iterations to reduce the diagonal elements to zero in this precision than in double precision.
        #     # an even more stringent requirement would be for the off diagonal elements also not to be changed by the rotation matrix but this can increase the number of iterations by 3-5 times for small 
        #     # matrices, maybe more for big matrices. In some extreme cases exiting before the off diagonal elements cease to be changed may result in more error in the singular vectors than desirable but I have not observed this yet.
        #     needToRotate = false
        #     return sinTheta,cosTheta,needToRotate
        # else
        #     needToRotate = true
        #     return sinTheta,cosTheta,needToRotate
        # end

        needToRotate = true
        return sinTheta, cosTheta, needToRotate
    else
        # off diagonal elements essentially zero
        sinTheta = 0
        cosTheta = 0
        needToRotate = false
        return sinTheta, cosTheta, needToRotate
    end
end

function oneSidedJacobiRotation!(A::AbstractArray{T,2}, j::Int, k::Int, sinTheta::T, cosTheta::T) where {T<:Real} # need AbstractArray instead of Array to handle Transpose, which is an AbstractArray but not an Array
    s = sinTheta
    c = cosTheta
    numrows, _ = size(A)

    for i in 1:numrows
        Aij = A[i, j]
        Aik = A[i, k]
        A[i, j] = c * Aij + s * Aik
        A[i, k] = -s * Aij + c * Aik
    end
end

"use this version if you want to handle exceeding maximum number of iterations. DESTRUCTIVE: matrix A is overwritten. powerPercentage must be be between 0 and 1"
function JacobiSolve(A::AbstractArray{T,2}, b::AbstractArray{T,1}; computeToLowRelativeAccuracy::Bool = false, conditionNumber::T = T(1e20))::Array{T,1} where {T<:Real}  # don't like having to declare this as Real rather than AbstractFloat because eps(Real) wouldn't necessarily be defined. But that's the way the IntervalArithmetic package is designed so to make it work with intervals have to this. The ForwardDiff package is also designed to take Real arguments. 
    local temp, total
    currentSum = zero(T)

    U, sigma, V, status = JacobiSVD(computeToLowRelativeAccuracy, A)
    temp = U' * b

    for i in 1:length(sigma)
        # println(S[0]/S[S.GetLength(0)-1]);
        if sigma[1] / sigma[i] > conditionNumber || sigma[i] == zero(T)
            temp[i] = zero(T)
        else
            temp[i] /= sigma[i]
        end
    end

    return V * temp

end

"DESTRUCTIVE: matrix A is overwritten. Returns U,sigma,V (not V'), iterationsexceeded   where iterationsexceeded will be either :exceededMax or :didNotExceedMax. If max iterations was exceeded it is possible the result is not accurate. But the results might also be fine."
function JacobiSVD(relativeaccuracy::Bool, A::AbstractArray{T,2}) where {T<:Real}
    acopy = copy(A)
    return JacobiSVD!(relativeaccuracy, acopy)
end

"sets flag to compute singular values to high relative accuracy. DESTRUCTIVE: matrix A is overwritten"
function JacobiSVD!(A::AbstractArray{T,2}) where {T<:Real}
    # out double[,] U, out double[] sigma, out double[,] V)
    return JacobiSVD!(false, A)
end

"DESTRUCTIVE: matrix A is overwritten. Returns U,sigma,V (not V'), iterationsexceeded   where iterationsexceeded will be either :exceededMax or :didNotExceedMax. If max iterations was exceeded it is possible the result is not accurate. But the results might also be fine."
function JacobiSVD!(computeToLowRelativeAccuracy::Bool, A::AbstractArray{T,2}) where {T<:Real}
    local status

    maxIterations = 35
    computingTranspose = false
    rowsA, colsA = size(A)

    # extend later
    if rowsA < colsA
        A = transpose(A)
        computingTranspose = true
        rowsA, colsA = size(A) # compute new size of transposed A
    end

    nrows = rowsA
    ncols = colsA
    n = colsA
    minRowsCols = min(nrows, ncols)
    sigma = Array{T,1}(undef, minRowsCols)
    localSigma = similar(sigma)

    numIterations = 0
    R = Matrix{T}(I, ncols, ncols)

    V = Array{T,2}(undef, n, n)
    U = Array{T,2}(undef, nrows, minRowsCols)
    while true

        rotatedThisIteration = false

        for j in 1:(n - 1)
            for k in (j + 1):n
                GTGjj = ATA(A, j, j)
                GTGjk = ATA(A, j, k)
                GTGkk = ATA(A, k, k)
                s, c, needToRotate = GivensRotation(computeToLowRelativeAccuracy, GTGjj, GTGjk, GTGkk)
                rotatedThisIteration = rotatedThisIteration || needToRotate
                if needToRotate
                    oneSidedJacobiRotation!(A, j, k, s, c)
                    oneSidedJacobiRotation!(R, j, k, s, c)
                end
            end
        end

        numIterations += 1
        if numIterations > maxIterations
            break
        end
        if !rotatedThisIteration
            break
        end
    end

    for i in 1:length(localSigma)
        localSigma[i] = 0
        for j in 1:rowsA
            # localSigma[i] += A[j, i] * A[j, i];
            localSigma[i] += A[j, i]^2 # use this instead of A[j,i]*A[j,i]. Interval arithmetic package will compute a much tighter inclusion
        end
        localSigma[i] = sqrt(localSigma[i])
    end

    # sort sigmas
    indices = sortperm(localSigma, rev = true)
    localSigma = localSigma[indices]


    for i in 1:minRowsCols
        sigma[i] = localSigma[i]  # copy just the true singular values

        for j in 1:rowsA
            if (sigma[i] != 0)
                U[j, i] = A[j, indices[i]] / sigma[i]
            else
                U[j, i] = 0
            end

        end
        for j in 1:n
            V[j, i] = R[j, indices[i]]
        end
    end

    if computingTranspose
        temp = V
        V = U
        U = temp
    end

    if numIterations > maxIterations
        status = :exceededMax
    else
        status = :didNotExceedMax
    end
    return (U = U, sigma = sigma, V = V, status = status)
end
