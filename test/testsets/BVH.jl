@testset "BVH" begin
    @testset "partition!" begin
        split = 0.5

        for i in 1:100000
            a = rand(5)
            b = copy(a)
            badresult::Bool = false

            lower, upper = partition!(a, (x) -> x, split)

            if lower !== nothing
                for i in lower
                    if i >= split
                        badresult = true
                        break
                    end
                end
            end

            if upper !== nothing
                for i in upper
                    if i <= split
                        badresult = true
                    end
                end
            end

            if badresult
                throw(ErrorException("array didn't partition: $(b)"))
            end

        end
    end
end # testset BVH
