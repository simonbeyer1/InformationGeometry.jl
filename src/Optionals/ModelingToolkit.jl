
using .ModelingToolkit

ExpTransform(Sys::ModelingToolkit.AbstractODESystem, idxs::AbstractVector{<:Bool}=trues(length(parameters(Sys))); kwargs...) = SystemTransform(Sys, exp, idxs; kwargs...)
LogTransform(Sys::ModelingToolkit.AbstractODESystem, idxs::AbstractVector{<:Bool}=trues(length(parameters(Sys))); kwargs...) = SystemTransform(Sys, log, idxs; kwargs...)
Power10Transform(Sys::ModelingToolkit.AbstractODESystem, idxs::AbstractVector{<:Bool}=trues(length(parameters(Sys))); kwargs...) = SystemTransform(Sys, x->(10.0)^x, idxs; kwargs...)
Log10Transform(Sys::ModelingToolkit.AbstractODESystem, idxs::AbstractVector{<:Bool}=trues(length(parameters(Sys))); kwargs...) = SystemTransform(Sys, log10, idxs; kwargs...)

"""
    SystemTransform(Sys::ODESystem, F::Function, idxs::AbstractVector{<:Bool}=trues(length(parameters(Sys)))) -> ODESystem
Transforms the parameters of a `ODESystem` according to `F`.
"""
function SystemTransform(Sys::ModelingToolkit.AbstractODESystem, F::Function, idxs::AbstractVector{<:Bool}=trues(length(parameters(Sys))))
    SubstDict = Dict(parameters(Sys) .=> [(idxs[i] ? F(x) : x) for (i,x) in enumerate(parameters(Sys))])
    NewEqs = [(equations(Sys)[i].lhs ~ substitute(equations(Sys)[i].rhs, SubstDict)) for i in 1:length(equations(Sys))]
    ODESystem(NewEqs, independent_variables(Sys)[1], states(Sys), parameters(Sys); name=nameof(Sys))
end



"""
    InformNames(DS::AbstractDataSet, sys::ODESystem, observables::AbstractVector{<:Int})
Copy the state names saved in `ODESystem` to `DS`.
"""
function InformNames(DS::AbstractDataSet, sys::ModelingToolkit.AbstractSystem, observables::Union{BoolArray,AbstractVector{<:Int}})
    newxnames = xnames(DS) == CreateSymbolNames(xdim(DS),"x") ? [string(ModelingToolkit.get_iv(sys))] : xnames(DS)
    newynames = ynames(DS) == CreateSymbolNames(ydim(DS),"y") ? string.(ModelingToolkit.get_states(sys)[observables]) : ynames(DS)
    InformNames(DS, newxnames, newynames)
end


# Use sys to infer state names of ODEsys
# Extend for other DEFunctions in the future
function DataModel(DS::AbstractDataSet, sys::ModelingToolkit.AbstractSystem, u0::Union{AbstractArray{<:Number},Function},
                        observables::Union{AbstractVector{<:Int},BoolArray,Function}=collect(1:length(u0)), args...; tol::Real=1e-7, Domain::Union{HyperCube,Nothing}=nothing, kwargs...)
    newDS = observables isa AbstractVector ? InformNames(DS, sys, observables) : DS
    DataModel(newDS, GetModel(sys, u0, observables; tol=tol, Domain=Domain, kwargs...), args...)
end


function GetModel(sys::ModelingToolkit.AbstractSystem, u0::Union{AbstractArray{<:Number},Function}, observables::Union{AbstractVector{<:Int},BoolArray,Function}=collect(1:length(u0));
                Domain::Union{HyperCube,Nothing}=nothing, inplace::Bool=true, kwargs...)
    # Is there some optimization that can be applied here? Modellingtoolkitize(sys) or something?
    # sys = Sys isa Catalyst.ReactionSystem ? convert(ODESystem, Sys) : Sys
    Model = if sys isa ModelingToolkit.AbstractODESystem
        GetModel(ODEFunction{inplace}(sys), u0, observables; Domain=Domain, inplace=inplace, kwargs...)
    else
        throw("Not programmed for $(typeof(sys)) yet, please convert to a ModelingToolkit.AbstractODESystem first.")
    end
    if Model isa ModelMap       Model = Model.Map    end
    pnames = ModelingToolkit.get_ps(sys) .|> string
    ylen = if observables isa Function      # ObservationFunction
        # Might still fail if states u are a Matrix.
        argnum = MaximalNumberOfArguments(observables)
        F = if argnum==1  z->observables(z)   elseif argnum==2  z->observables(z,0.1)
            elseif argnum==3 z->observables(z,0.1,GetStartP(length(pnames)))    else throw("Error") end
        num = GetArgLength(F)
        length(F(ones(num)))
    else
        observables isa BoolArray ? sum(observables) : length(observables)
    end
    plen = if Domain isa HyperCube
        length(Domain)
    elseif u0 isa AbstractArray     # Vector / Matrix
        # initial conditions given as array means the parameters are only the ps in sys
        length(pnames)
    else        # SplitterFunction
        # May well fail depending on how splitter function is implemented
        GetArgLength(u0)
    end
    xyp = (1, ylen, plen)
    Domain = isa(Domain, Bool) ? FullDomain(xyp[3], 1e5) : Domain

    pnames = plen - length(pnames) > 0 ? vcat(CreateSymbolNames(plen - length(pnames), "u"), pnames) : pnames
    # new(Map, InDomain, Domain, xyp, pnames, StaticOutput, inplace, CustomEmbedding)
    ModelMap(Model, nothing, Domain, xyp, pnames, Val(false), Val(false), Val(true))
end

@info "InformationGeometry: Loaded ModelingToolkit extensions."
