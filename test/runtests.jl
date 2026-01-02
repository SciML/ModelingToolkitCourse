using Test

@testset "ModelingToolkitCourse.jl" begin
    @testset "Explicit Imports" begin
        include("explicit_imports.jl")
    end
end
