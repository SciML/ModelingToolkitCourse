using ExplicitImports
using ModelingToolkitCourse
using Test

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(ModelingToolkitCourse) === nothing
    @test check_no_stale_explicit_imports(ModelingToolkitCourse) === nothing
end
