using ModelingToolkit, OrdinaryDiffEq, Test

@parameters t σ ρ β
@variables x(t) y(t) z(t) F(t) u(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x) + F,
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

aliases = [u ~ x + y - z]
lorenz1 = ODESystem(eqs,inputs=[F],outputs=aliases,name=:lorenz1)
lorenz2 = ODESystem(eqs,inputs=[F],outputs=aliases,name=:lorenz2)

connections = [lorenz1.F ~ lorenz2.u,
               lorenz2.F ~ lorenz1.u]
connected = ODESystem(Equation[],t,[],[],outputs=connections,systems=[lorenz1,lorenz2])

sys = connected

@variables lorenz1₊F lorenz2₊F
@test inputs(connected) == Variable[lorenz1₊F, lorenz2₊F]
@test isequal(outputs(connected),
              [connections...,
               lorenz1.u ~ lorenz1.x + lorenz1.y - lorenz1.z,
               lorenz2.u ~ lorenz2.x + lorenz2.y - lorenz2.z])

collapsed_eqs = [D(lorenz1.x) ~ (lorenz1.σ * (lorenz1.y - lorenz1.x) +
                                 (lorenz2.x + lorenz2.y - lorenz2.z)),
                 D(lorenz1.y) ~ lorenz1.x * (lorenz1.ρ - lorenz1.z) - lorenz1.y,
                 D(lorenz1.z) ~ lorenz1.x * lorenz1.y - (lorenz1.β * lorenz1.z),
                 D(lorenz2.x) ~ (lorenz2.σ * (lorenz2.y - lorenz2.x) +
                                 (lorenz1.x + lorenz1.y - lorenz1.z)),
                 D(lorenz2.y) ~ lorenz2.x * (lorenz2.ρ - lorenz2.z) - lorenz2.y,
                 D(lorenz2.z) ~ lorenz2.x * lorenz2.y - (lorenz2.β * lorenz2.z)]

simplifyeqs(eqs) = Equation.((x->x.lhs).(eqs), simplify.((x->x.rhs).(eqs)))

@test isequal(simplifyeqs(equations(connected)), simplifyeqs(collapsed_eqs))
