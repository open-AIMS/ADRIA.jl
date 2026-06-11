using Pkg
Pkg.activate(".")
try
    using ADRIA
catch e
    showerror(stdout, e, catch_backtrace())
end
