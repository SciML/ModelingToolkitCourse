using Documenter

pages = [
    "Home" => "index.md",
    "lectures/lecture1.md",
    "lectures/lecture2.md",
]

makedocs(sitename = "ModelingToolkit Course",
    authors = "Chris Rackauckas",
    modules = [],
    clean = true, doctest = false, linkcheck = true,
    warnonly = [:missing_docs],
    format = Documenter.HTML(assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/ModelingToolkitCourse/stable/"),
    pages = pages)

deploydocs(repo = "github.com/SciML/ModelingToolkitCourse.git";
    push_preview = true)
