function rhosp_mixed(probs::Array{Float64,1}, states::Union{Array{State,1},Array{State_fixed,1}})
    if round(sum(probs), digits=15) != 1
        throw(ArgumentError("Your probabilities must add up to 1!"))
    end

    ls = length(states)
    if length(probs) != ls
        throw(ArgumentError("There must be the same ammount of probabilities and states!"))
    end

    #for i in 1:length(states)
    #    stat = st(states[i])
    #    if round(stat'*stat, digits=15) != 1
    #        throw(ArgumentError("All states must be normalized to 1!"))
    #    end
    #end

    rhos = [probs[i]*rhosp(states[i]) for i in 1:ls]
    rho = sum(rhos)
end


function eigensp_mixed(probs::Array{Float64,1}, states::Union{Array{State,1},Array{State_fixed,1}})
    if round(sum(probs), digits=15) != 1
        throw(ArgumentError("Your probabilities must add up to 1!"))
    end

    ls = length(states)
    if length(probs) != ls
        throw(ArgumentError("There must be the same ammount of probabilities and states!"))
    end

    #for i in 1:length(states)
    #    stat = st(states[i])
    #    if round(stat'*stat, digits=15) != 1
    #        throw(ArgumentError("All states must be normalized to 1!"))
    #    end
    #end

    eig = sort(eigvals(Matrix(rhosp_mixed(probs, states))), rev=true)
    return eig
end
