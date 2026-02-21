"""
    generate_ODENLStepData(sys, u0, p, mm, nlstep_compile, nlstep_scc)

Generate the NLStep data for implicit ODE solvers. This is a stub that throws an error
if called without ModelingToolkit loaded. The actual implementation is provided by
ModelingToolkit when it is loaded.
"""
function generate_ODENLStepData(sys, u0, p, mm, nlstep_compile, nlstep_scc)
    error(
        """
        `nlstep=true` requires ModelingToolkit.jl to be loaded.
        Please add `using ModelingToolkit` to your code before creating an ODEProblem with `nlstep=true`.
        """
    )
end

Base.@nospecializeinfer @fallback_iip_specialize function SciMLBase.ODEFunction{iip, spec}(
        sys::System; @nospecialize(u0 = nothing), @nospecialize(p = nothing), tgrad = false, jac = false,
        t = nothing, eval_expression = false, eval_module = @__MODULE__, sparse = false,
        steady_state = false, checkbounds = false, sparsity = false, @nospecialize(analytic = nothing),
        simplify = false, cse = true, @nospecialize(initialization_data = nothing), expression = Val{false},
        check_compatibility = true, nlstep = false, nlstep_compile = true, nlstep_scc = false,
        optimize = nothing, kwargs...
    ) where {iip, spec}
    check_complete(sys, ODEFunction)
    check_compatibility && check_compatible_system(ODEFunction, sys)

    f = generate_rhs(
        sys; expression, wrap_gfw = Val{true},
        eval_expression, eval_module, checkbounds = checkbounds, cse,
        optimize, kwargs...
    )

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on ODEFunction.")
        end
        if expression == Val{true}
            f = :($(SciMLBase.wrapfun_iip)($f, ($u0, $u0, $p, $t)))
        else
            f = SciMLBase.wrapfun_iip(f, (u0, u0, p, t))
        end
    end

    if tgrad
        _tgrad = generate_tgrad(
            sys; expression, wrap_gfw = Val{true},
            simplify, cse, eval_expression, eval_module, checkbounds, optimize, kwargs...
        )
    else
        _tgrad = nothing
    end

    if jac
        _jac = generate_jacobian(
            sys; expression, wrap_gfw = Val{true},
            simplify, sparse, cse, eval_expression, eval_module, checkbounds, optimize,
            kwargs...
        )
    else
        _jac = nothing
    end

    M = calculate_massmatrix(sys)
    _M = concrete_massmatrix(M; sparse, u0)

    if nlstep
        ode_nlstep = generate_ODENLStepData(sys, u0, p, M, nlstep_compile, nlstep_scc)
    else
        ode_nlstep = nothing
    end

    observedfun = ObservedFunctionCache(
        sys; expression, steady_state, eval_expression, eval_module, checkbounds, cse, optimize
    )

    _W_sparsity = W_sparsity(sys)
    W_prototype = calculate_W_prototype(_W_sparsity; u0, sparse)

    args = (; f)
    kwargs = (;
        sys = sys,
        jac = _jac,
        tgrad = _tgrad,
        mass_matrix = _M,
        jac_prototype = W_prototype,
        observed = observedfun,
        sparsity = sparsity ? _W_sparsity : nothing,
        analytic = analytic,
        initialization_data,
        nlstep_data = ode_nlstep,
    )

    odefn = maybe_codegen_scimlfn(expression, ODEFunction{iip, spec}, args; kwargs...)
    # Erase the OverrideInitData and ODENLStepData types so that all
    # AutoSpecialize ODEFunctions share the same type regardless of
    # model-specific init functions.
    if expression != Val{true} && spec === SciMLBase.AutoSpecialize
        odefn = _erase_init_data_type(odefn)
    end
    return odefn
end

"""
Reconstruct the ODEFunction with abstract union types for the `initialization_data` and
`nlstep_data` type parameters. This ensures all AutoSpecialize ODEFunctions have identical
types, preventing recompilation of `promote_f` and solver code for each model.
"""
function _erase_init_data_type(f::SciMLBase.ODEFunction)
    return SciMLBase.ODEFunction{
        SciMLBase.isinplace(f), SciMLBase.specialization(f),
        typeof(f.f), typeof(f.mass_matrix),
        typeof(f.analytic), typeof(f.tgrad),
        typeof(f.jac), typeof(f.jvp), typeof(f.vjp), typeof(f.jac_prototype),
        typeof(f.sparsity), typeof(f.Wfact), typeof(f.Wfact_t), typeof(f.W_prototype),
        typeof(f.paramjac),
        typeof(f.observed), typeof(f.colorvec),
        typeof(f.sys),
        Union{Nothing, SciMLBase.OverrideInitData},
        Union{Nothing, SciMLBase.ODENLStepData},
    }(
        f.f, f.mass_matrix, f.analytic, f.tgrad, f.jac,
        f.jvp, f.vjp, f.jac_prototype, f.sparsity, f.Wfact,
        f.Wfact_t, f.W_prototype, f.paramjac,
        f.observed, f.colorvec, f.sys, f.initialization_data, f.nlstep_data
    )
end

Base.@nospecializeinfer @fallback_iip_specialize function SciMLBase.ODEProblem{iip, spec}(
        sys::System, @nospecialize(op), tspan;
        @nospecialize(callback = nothing), check_length = true, eval_expression = false,
        expression = Val{false}, eval_module = @__MODULE__, check_compatibility = true,
        _skip_events = false, kwargs...
    ) where {iip, spec}
    check_complete(sys, ODEProblem)
    check_compatibility && check_compatible_system(ODEProblem, sys)

    f, u0,
        p = process_SciMLProblem(
        ODEFunction{iip, spec}, sys, op;
        t = tspan !== nothing ? tspan[1] : tspan, check_length, eval_expression,
        eval_module, expression, check_compatibility, kwargs...
    )

    kwargs = process_kwargs(
        sys; expression, callback, eval_expression, eval_module, op, _skip_events, kwargs...
    )

    ptype = getmetadata(sys, ProblemTypeCtx, StandardODEProblem())
    args = (; f, u0, tspan, p, ptype)
    maybe_codegen_scimlproblem(expression, ODEProblem{iip}, args; kwargs...)
end

@fallback_iip_specialize function DiffEqBase.SteadyStateProblem{iip, spec}(
        sys::System, op; check_length = true, check_compatibility = true,
        expression = Val{false}, kwargs...
    ) where {iip, spec}
    check_complete(sys, SteadyStateProblem)
    check_compatibility && check_compatible_system(SteadyStateProblem, sys)

    f, u0,
        p = process_SciMLProblem(
        ODEFunction{iip}, sys, op;
        steady_state = true, check_length, check_compatibility, expression,
        is_steadystateprob = true, kwargs...
    )

    kwargs = process_kwargs(sys; expression, kwargs...)
    args = (; f, u0, p)

    maybe_codegen_scimlproblem(expression, SteadyStateProblem{iip}, args; kwargs...)
end

function check_compatible_system(
        T::Union{
            Type{ODEFunction}, Type{ODEProblem}, Type{DAEFunction},
            Type{DAEProblem}, Type{SteadyStateProblem},
        },
        sys::System
    )
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    return check_is_continuous(sys, T)
end
