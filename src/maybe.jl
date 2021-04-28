# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

export Maybe, MaybeVector

"""
Haskell-style [option type](https://en.wikipedia.org/wiki/Option_type) primarily used to allow optional inputs/outputs
for functions.

For example, in `src/Examples/docs_examples.jl`, we have `draw_cooketriplet(filename::Maybe{AbstractString} = nothing)`.
If the function is called with a `filename` argument, then it is passed as a string, and `Vis.save(filename)` saves an
image of the drawn example. If there is no argument, then the save call is dispatched to the empty `function
Vis.save(::Nothing) end`. This allows the user to run the example without saving a file.
"""
const Maybe{T} = Union{T, Nothing}

"""
Extension of `Maybe{T}` which allows things like `[nothing, 1, 1.0, nothing] isa MaybeVector{Real}`
"""
const MaybeVector{T} = Vector{<:Maybe{T}}
