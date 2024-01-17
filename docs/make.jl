using Documenter, ModelingToolkitCourse

pages = [
    "Home" => "index.md",
    "Syllabus" => "syllabus.md",
    "lectures/lecture1.md",
    "lectures/lecture2.md",
    "lectures/lecture3.md",
    "lectures/lecture6.md",
    "lectures/lecture7.md",
]

ENV["GKSwstype"] = "100"
using Plots

makedocs(sitename = "ModelingToolkit Course",
    authors = "Chris Rackauckas",
    modules = [ModelingToolkitCourse],
    clean = true, doctest = false, linkcheck = true,
    warnonly = [:missing_docs],
    format = Documenter.HTML(assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/ModelingToolkitCourse/stable/"),
    pages = pages)

deploydocs(repo = "github.com/SciML/ModelingToolkitCourse.git";
    push_preview = true)
