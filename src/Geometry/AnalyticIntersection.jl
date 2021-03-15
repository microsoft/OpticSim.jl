# These constants determine the order in which the g coefficients are stored in the nullspace matrix, and the corresponding moving line matrix. Do not change these. Every other function depends on this.
# later will change code so changing these constants will automagically work.
const lineindex = 3
const spatialindex = 1
const orderindex = 2

"spatial dimension of curve represented as an array of coefficients `x[i] = ∑Bj(θ)*x[i,j]` where `Bj(θ)` is the curve basis"
curvedimension(x::Array) = size(x, spatialindex)

"highest polynomial power of the curve represented as an array of coefficients `x[i] = ∑Bj(θ)*x[i,j]` where `Bj(θ)` is the curve basis"
curveorder(x::Array) = size(x, orderindex) - 1

numcoefficients(x::Array) = size(x, orderindex)

"spatial dimension of the moving line represented as an array of coefficients `g[i] = ∑Bl(θ)*gl[i,j]` where `Bl(θ)` is the polynomial basis"
linedimension(x::Array) = size(x, 1) - 1

"number of lines in moving line array"
numberoflines(x::Array) = size(x, lineindex)

movingline(x::Array, linenum::Integer) = x[:, :, linenum]

"returns a matrix expressing the relationship `[x(θ) 1]⋅g(θ) = 0`. The vectors in the right nullspace of this matrix contain the coefficients of the moving lines `gᵢ(θ)`."
function orthogonalitymatrix(curve::Array{T,2}, movinglineorder) where {T}
    x = curve
    qx = curveorder(x)
    qg = movinglineorder
    d = curvedimension(x)

    C = zeros(T, qx + qg + 1, (d + 1) * (qg + 1))

    for i in 1:d, j in 1:(qx + 1), l in 1:(qg + 1)
        C[(j - 1) + (l - 1) + 1, (i - 1) * (qg + 1) + l] = x[i, j]
    end

    # now set the final term which does not include any components of the curve x(theta)
    for row in 1:(qg + 1)
        col = d * (qg + 1) + row
        C[row, col] = one(T)
    end

    return C
end

"returns 3D array indexed like this: `x[line curve order,spatial dimension, line number]``"
function extractmovinglines(Vt, nullspacesize, movinglineorder, dimension)
    rows, _ = size(Vt)
    # temp = Vt[(rows-nullspacesize + 1):end,:]
    # println("rows $rows null $nullspacesize")
    # show(IOContext(stdout), "text/plain", temp)

    lines = permutedims(Vt[(rows - nullspacesize + 1):end, :]) # use this instead of transpose because transpose is recursive so it tries to transpose the elements of the array and crashes for Expr element types.
    # show(IOContext(stdout), "text/plain", lines)
    linecoeffs = movinglineorder + 1
    return permutedims(reshape(lines, linecoeffs, (dimension + 1), :), (2, 1, 3))
end

"`movinglines[:,i]` is the ith moving line. For `li = movinglines[:,i] (dimension+1,lineorder) = size(li)`. `rline[:,1] = pt1` and `rline[:,2] = pt2`. The line equation is `pt1 + alpha*pt2`."
function matricesforeigen(ray::AbstractRay{T,N}, movinglines::Array{T,3}) where {T,N}
    pt1 = origin(ray)
    pt2 = direction(ray)

    dimension = size(pt1, 1)
    numcoeffs = numcoefficients(movinglines)
    numlines = numberoflines(movinglines)

    A = Array{T,2}(undef, numcoeffs, numlines)
    B = Array{T,2}(undef, numcoeffs, numlines)

    for i in 1:numlines
        li = movinglines[:, :, i]

        for l in 1:numcoeffs
            asum = zero(T)
            bsum = zero(T)

            for j in 1:dimension
                asum += pt1[j] * li[j, l]
                bsum += pt2[j] * li[j, l]
            end

            asum += li[dimension + 1, l] # add in the final term of the sum that is not dotted with rline

            bsum = -bsum # reverse sign to get matrices in A-alpha*B form
            A[l, i] = asum
            B[l, i] = bsum
        end
    end
    return (A, B)
end

function matrixsizes(dimension, movinglineorder, curveorder)
    cols, rows = ((dimension + 1) * (movinglineorder + 1), curveorder + movinglineorder + 1)  # cols is the number of coefficients for the moving line. The line has d+1 terms for the line(plane/hyperplane) equation and movinglineorder+1 coefficients for the curve describing the moving line.
    nullspacesize = cols - rows

    return cols, rows, nullspacesize
end

"Evaluates a curve defined in the power basis. Curves and moving lines accessed like this: `[xi,ci]` where `xi` is the dimension index, and `ci` is the coefficient index."
function evaluatecurve(x::Array{T,2}, theta::Real) where {T<:Real}
    # dim,coefficients = size(x)
    dim = curvedimension(x)
    numcoeffs = numcoefficients(x)
    result = Array{T,1}(undef, dim)

    for dimension in 1:dim
        power = one(T)
        sum = zero(T)
        for coefficient in 1:numcoeffs
            sum += x[dimension, coefficient] * power
            power *= theta
        end
        result[dimension] = sum
    end
    return result
end

function validintersection(rline::AbstractRay{T,N}, linealpha, curve::Array{T,2}, curvetheta) where {T<:Real,N}
    if linealpha < 0
        return false
    end
    linepoint = point(rline, linealpha)
    curvepoint = evaluatecurve(curve, curvetheta)
    return isapprox(linepoint, curvepoint, rtol = 1e-10)
end

function eigenresults(rline::AbstractRay{T,N}, curve::Array{T,2}) where {T<:Real,N} # force rline and curve to use same number type to avoid expensive runtime conversion.
    dim = curvedimension(curve)
    orderofcurve = curveorder(curve)
    movinglineorder = max(orderofcurve, Int64(ceil(orderofcurve * 2 / dim - 1))) # Hve two constraints for movingline order. Must be big enough to result in a nullspace >= orderofcurve for the first step of computing the moving lines. Also need need at least as many eigenvalues as possible intersections of the line with the curve. Eigenmatrix is square of size qgxqg so need this matrix to be at least of size qx by qx. Take the greater of the the two constraint values.

    cmat = orthogonalitymatrix(curve, movinglineorder)

    fact = svd(cmat, full = true)
    nullspace = fact.Vt

    _, _, nullspacesize = matrixsizes(dim, movinglineorder, orderofcurve)

    movinglines = extractmovinglines(nullspace, nullspacesize, movinglineorder, dim)
    A, B = matricesforeigen(rline, movinglines)

    eig = eigen(A[:, 1:4]', B[:, 1:4]') # transpose so [A-lambdaB]'*Pl = 0

    return eig
end

"returns an array of intersection points. Each element in the array is (`[x,y,...],alpha,theta)` where `[x,y,...]` is the n-dimensional intersection point, alpha is the line parameter value at the intersection point, and theta is the curve parameter value at the intersection point"
function intersections(rline::AbstractRay{T,N}, curve::Array{T,2})::Array{Tuple{Array{T,1},T,T}} where {T<:Real,N}
    eigenstuff = eigenresults(rline, curve)
    result = Array{Tuple{Array{T,1},T,T}}(undef, 0)

    for index in CartesianIndices(eigenstuff.values)
        eigenvector = view(eigenstuff.vectors, :, index)
        eigenvalue = eigenstuff.values[index]
        # only take the real part of the solution. Due to roundoff might have real roots with very small complex part.
        theta = real(eigenvector[2] / eigenvector[1]) # this is only correct for curves defined in power basis. Once new Curve types are defined this function will be determined by the curve basis of the curve. Something like curveParameter(::CurveBasis,eigenvector)
        println(theta)
        if validintersection(rline, eigenvalue, curve, theta)
            push!(result, (point(rline, eigenvalue), eigenvalue, theta)) # return point on curve and alpha and theta parametric values corresponding to this point.
        end
    end

    return result
end
export intersections
