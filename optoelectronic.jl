# <md>
# chimera states in a small ring of electro-optic oscillators

# In[]
# loading necessary packages
#using DifferentialEquations
using Plots; pyplot()
using LightGraphs
#using GraphPlot
using LinearAlgebra
using Statistics
using BenchmarkTools
using ProgressMeter
using LaTeXStrings
using DelimitedFiles

# In[]
n = 6   # size of the ring
k = 2    # number of neighbors coupled (in each direction)
σ = -.7    # coupling strength
β = 2*π/3 - 4*σ   # parameter of the electro-optic oscillators
ξ = π/6     # offset parameter
r = 0     # rewiring probability
ring = watts_strogatz(n, 2*k, r)
adj = Matrix(adjacency_matrix(ring))

# In[]
# mapping function for the electro-optic oscillators
function I(x)
    return (1-cos(x))/2
end

# In[]
# simulating the full system to study chimera states
T = 1000000 # total simulation time
x = zeros(n,T) # oscillator states
x[:,1] = 2*π*[.1,.5,.5,.7,.5,.5] + 1e-3*rand(n,1) # initial condition close to a chimera state
p = Progress(T, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow) # progress bar
@time for t in 2:T
    for i in 1:n
        x[i,t] = β*I(x[i,t-1]) + σ*adj[i,:]⋅I.(x[:,t-1]) + ξ
    end
    x[:,t] = mod.(x[:,t],2*π)
    next!(p)
end


# In[]
# plot showing the states of the oscillators in a chimera state (coherent cluster: green; incoherent cluster: orange and magenta)
t = 500 # show the final t iterations
scatter(T-t:T,x[1,T-t:T],marker=(:circle,3,.6,:magenta,stroke(0,.2,:black)),background_color=RGB(0.2, 0.2, 0.2),legend=false)
scatter!(T-t:T,x[4,T-t:T],marker=(:circle,3,.6,:orange,stroke(0,.2,:black)))
for i = [2,3,5,6]
    scatter!(T-t:T,x[i,T-t:T],marker=(:circle,3,.6,stroke(1,.2,:white)))
end
title!("Time Series for Chimeras")
xaxis!(L"t",font(15))
yaxis!(L"x_i",font(15),(-.1,2*π))
yticks!([0,π,2*π],[L"0",L"\pi",L"2\pi"])


# In[]
# plot showing the input from the incoherent cluster to the coherent cluster
@show mean(σ*I.(x[1,:])+σ*I.(x[4,:])) # mean of the input
@show std(σ*I.(x[1,:])+σ*I.(x[4,:])) # s.t.d. of the input
t = 500 # show the final t iterations
scatter(T-t:T,σ*I.(x[1,T-t:T])+σ*I.(x[4,T-t:T]),marker=(:circle,3,.6,:red,stroke(0,.2,:black)),background_color=RGB(0.2, 0.2, 0.2),legend=false)
title!("Time Series of the Driving")
xaxis!(L"t",font(15))
yaxis!("Input",font(15))


# In[]
# simulating the two-oscillator reduced system
T = 1000000
w = zeros(2,T) # oscillator states
w[:,1] = 2*π*rand(2,1) # random initial condition
p = Progress(T, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
@time for t in 2:T
    # quotient of the 4-node cluster
    w[1,t] = β*I(w[1,t-1]) + 2*σ*(I(w[2,t-1])) + 2*σ*(I(w[1,t-1]))  + ξ
    # quotient of the 2-node cluster
    w[2,t] = β*I(w[2,t-1]) + 4*σ*(I(w[1,t-1])) + ξ
    w[:,t] = mod.(w[:,t],2*π)
    next!(p)
end

# In[]
# plot showing the states of the oscillators in the two-oscillator reduced system (which is unstable in the full system)
t = 500 # show the final t iterations
scatter(T-t:T,w[1,T-t:T],marker=(:circle,3,.6,:orange,stroke(0,.2,:black)),background_color=RGB(0.2, 0.2, 0.2),legend=false)
scatter!(T-t:T,w[2,T-t:T],marker=(:circle,3,.6,:green,stroke(0,.2,:black)),background_color=RGB(0.2, 0.2, 0.2),legend=false)
title!("Reduced System")
xaxis!(L"t",font(15))
yaxis!(L"x_i",font(15),(-.1,2*π))
yticks!([0,π,2*π],[L"0",L"\pi",L"2\pi"])

# In[]
# plot showing the input from one cluster to the other
@show mean(2*σ*I.(w[2,:])) # mean of the input
@show std(2*σ*I.(w[2,:])) # s.t.d. of the input
t = 500 # show the final t iterations
scatter(T-t:T,2*σ*I.(w[2,T-t:T]),marker=(:circle,3,.6,:red,stroke(0,.2,:black)),background_color=RGB(0.2, 0.2, 0.2),legend=false)
title!("Time Series of the Driving")
xaxis!(L"t",font(15))
yaxis!("Input",font(15))
