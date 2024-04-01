using Documenter, ModelingToolkitCourse

pages = [
    "Home" => "index.md",
    "lectures/lecture1.md",
    "lectures/lecture2.md",
    "lectures/lecture3.md",
    "lectures/lecture4.md",
    "lectures/lecture6.md",
    "lectures/lecture7.md",
    "lectures/lecture8.md",
]

ENV["GKSwstype"] = "100"
using Plots

makedocs(sitename = "ModelingToolkit Course",
    authors = "Chris Rackauckas",
    modules = [ModelingToolkitCourse],
    clean = true, doctest = false, linkcheck = true,
    linkcheck_ignore = ["https://epubs.siam.org/doi/10.1137/0903023",
    "https://link.springer.com/book/10.1007/978-3-642-05221-7",
    "http://www.siam.org/journals/auth-info.php"],
    format = Documenter.HTML(assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/ModelingToolkitCourse/stable/"),
    pages = pages)

deploydocs(repo = "github.com/SciML/ModelingToolkitCourse.git";
    push_preview = true)
