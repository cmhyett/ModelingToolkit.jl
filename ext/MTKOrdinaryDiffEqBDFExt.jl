module MTKOrdinaryDiffEqBDFExt

using ModelingToolkit
using OrdinaryDiffEqBDF
using PrecompileTools
using ModelingToolkit: t_nounits, D_nounits

@setup_workload begin
    @variables x(t_nounits)
    prob = ODEProblem(
        mtkcompile(System([D_nounits(x) ~ 2x + 1], t_nounits; name = :precompile_bdf)),
        [x => 1.0], (0.0, 1.0))
    @compile_workload begin
        solve(prob, FBDF())
    end
end

end
