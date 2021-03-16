struct VariableUnit end
struct VariableConnectType end
Symbolics.option_to_metadata_type(::Val{:unit}) = VariableUnit
Symbolics.option_to_metadata_type(::Val{:connect}) = VariableConnectType

"""
$(SIGNATURES)

Takes a list of pairs of `variables=>values` and an ordered list of variables
and creates the array of values in the correct order with default values when
applicable.
"""
function varmap_to_vars(varmap, varlist; defaults=Dict())
    # Edge cases where one of the arguments is effectively empty.
    if varmap isa DiffEqBase.NullParameters || varmap === nothing || isempty(varmap)
        if isempty(defaults)
            isempty(varlist) || throw_missingvars(varlist)
            return nothing
        else
            varmap = Dict()
        end
    end

    T = typeof(varmap)
    # We respect the input type
    container_type = T <: Dict ? Array : T

    if eltype(varmap) <: Pair # `varmap` is a dict or an array of pairs
        varmap = todict(varmap)
        rules = Dict(varmap)
        vals = _varmap_to_vars(varmap, varlist; defaults=defaults)
    else # plain array-like initialization
        vals = varmap
    end

    if isempty(vals)
        return nothing
    elseif container_type <: Tuple
        (vals...,)
    else
        SymbolicUtils.Code.create_array(container_type, eltype(vals), Val(length(vals)), vals...)
    end
end

function _varmap_to_vars(varmap::Dict, varlist; defaults=Dict())
    varmap = merge(defaults, varmap) # prefers the `varmap`
    varmap = Dict(Symbolics.diff2term(value(k))=>value(varmap[k]) for k in keys(varmap))
    # resolve symbolic parameter expressions
    for (p, v) in pairs(varmap)
        varmap[p] = fixpoint_sub(v, varmap)
    end
    T′ = eltype(values(varmap))
    T = Base.isconcretetype(T′) ? T′ : Base.promote_typeof(values(varmap)...)
    out = Vector{T}(undef, length(varlist))
    missingvars = setdiff(varlist, keys(varmap))
    isempty(missingvars) || throw_missingvars(missingvars)

    for (i, var) in enumerate(varlist)
        out[i] = varmap[var]
    end
    out
end

@noinline throw_missingvars(vars) = throw(ArgumentError("$vars are missing from the variable map."))
