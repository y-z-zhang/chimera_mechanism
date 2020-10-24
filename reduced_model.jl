# <md>
# strong chimeras in two-cluster networks

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
n = 3     # size of the clusters
c = .2    # intra- and inter-coupling strength ratio
r = 1.2   # parameter of the optoelectronic map
K = 1.2   # overall coupling strength

# In[]
# discrete map model of optoelectronic oscillators
function f(x)
    return sin(x+π/4)^2
end

# In[]
# reduce the two-cluster system to a two-oscillator system
# simulate the two-oscillator reduced system to study the two-cluster coherent state
T = 100000 # total simulation time
x = zeros(2,T) # oscillator states
x[:,1] = rand(2,1) # random initial condition
p = Progress(T, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow) # progress bar
@time for t in 2:T
    x[1,t] = r*f(x[1,t-1]) + c*n*K*(f(x[2,t-1])-f(x[1,t-1]))
    x[2,t] = r*f(x[2,t-1]) + c*n*K*(f(x[1,t-1])-f(x[2,t-1]))
    next!(p)
end


# In[]
# plot showing the states of the two oscillators in the reduced system
scatter(1:T,x[1,:],marker=(:circle,4,.6,:blue,stroke(0,.2,:black)),background_color=RGB(0.2, 0.2, 0.2),legend=false)
scatter!(1:T,x[2,:],marker=(:circle,3,.6,:green,stroke(0,.2,:black)),background_color=RGB(0.2, 0.2, 0.2),legend=false)
title!("Reduced System for Two-Cluster Coherent States")
xaxis!(L"t",font(15),(T-99,T))
yaxis!(L"x_i",font(15),(1,1.1))
#yticks!([-π,0,π],[L"\pi",L"0",L"\pi"])


# In[]
# plot showing the effective input from one cluster to the other
scatter(1:T,c*n*K*(f.(x[2,:])-f.(x[1,:])),marker=(:circle,3,.6,:red,stroke(0,.2,:black)),background_color=RGB(0.2, 0.2, 0.2),legend=false)
title!("Reduced System for Two-Cluster Coherent States")
xaxis!(L"t",font(15),(T-99,T))
yaxis!("Input",font(15),(-.05,.05))


# In[]
# create the quotient network that corresponds to chimera states
function quotient_net(c)
    ring = watts_strogatz(n, 2, 0)
    adj = Matrix(adjacency_matrix(ring))
    return adj = [adj c*n*ones(n,1); c*ones(1,n) 0]
end


# In[]
# simulate the (n+1)-oscillator reduced system to study chimera states
function chimera(adj,T,IC)
    x = zeros(n+1,T) # oscillator states
    driving = zeros(T) # effective input from the incoherent cluster to the coherent cluster
    x[:,1] = IC # initial condition
    for t in 2:T
        for i in 1:n # updating oscillators in the incoherent cluster
            x[i,t] = r*f(x[i,t-1]) + K*adj[i,:]⋅(f.(x[:,t-1]).-f(x[i,t-1]))
        end
        driving[t] = K*adj[n+1,:]⋅(f.(x[:,t-1]).-f(x[n+1,t-1]))
        x[n+1,t] = r*f(x[n+1,t-1]) + driving[t] # updating the oscillator representing the coherent cluster
    end
    # return the mean and the s.t.d. of the effective input from the incoherent cluster to the coherent cluster, and the final oscillator states
    return mean(driving[T-10000:T]), std(driving[T-10000:T]), x[:,T]
end

# In[]
# calculate the mean and the s.t.d. of the effective input for different values of c
m = 100 # number of different values of c
T = 100000 # total simulation time for each c
μ = zeros(m) # mean of the effective input
σ = zeros(m) # s.t.d. of the effective input
IC = r*rand(n+1) # random initial conditions
c_range = range(0,stop=.4,length=m)
p = Progress(m, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
for i = 1:m
    c = c_range[i]
    adj = quotient_net(c)
    μ[i], σ[i], IC1 = chimera(adj,T,IC)
    next!(p)
end

# In[]
writedlm("input.txt",[μ,σ],',')

# In[]
# plot showing the mean and the s.t.d. of the effective input for different values of c
scatter(c_range,μ,marker=(:circle,3,.6,:green,stroke(0,.2,:black)),background_color=RGB(0.2, 0.2, 0.2),label = "μ")
scatter!(c_range,σ,marker=(:circle,3,.6,:blue,stroke(0,.2,:black)),background_color=RGB(0.2, 0.2, 0.2),label = "σ")
title!("Reduced System for Chimera States")
xaxis!(L"c",font(15))
yaxis!(L"\mu,\sigma",font(15))

# In[]
# simulate the (n+1)-oscillator reduced system to obtain a typical time series of the effective input
c = .2
adj = quotient_net(c)
T = 1000000
x = zeros(n+1,T) # oscillator states
driving = zeros(T) # effective input
x[:,1] = rand(n+1,1) # random initial condition
p = Progress(T, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
for t in 2:T
    for i in 1:n
        x[i,t] = r*f(x[i,t-1]) + K*adj[i,:]⋅(f.(x[:,t-1]).-f(x[i,t-1]))
    end
    driving[t] = K*adj[n+1,:]⋅(f.(x[:,t-1]).-f(x[n+1,t-1]))
    x[n+1,t] = r*f(x[n+1,t-1]) + driving[t]
    next!(p)
end
writedlm("driving.txt",driving,',') # time series of the effective input
@show mean(driving[T-100000:T]), std(driving[T-100000:T])

# In[]
# plot showing states of the (n+1) oscillators in the reduced system (coherent cluster: blue; incoherent cluster: green, red, magenta)
scatter(T-199:T,x[n+1,T-199:T],marker=(:circle,3,.6,:blue,stroke(1,.2,:white)),background_color=RGB(0.2, 0.2, 0.2),legend=false)
for i = 1:n
    scatter!(T-199:T,x[i,T-199:T],marker=(:circle,3,.6,stroke(0,.2,:black)))
end
title!("Reduced System for Chimera States")
xaxis!(L"t",font(15))
yaxis!(L"x_i",font(15))

# In[]
# plot showing a typical time series of the effective input
scatter(T-1000:T,driving[T-1000:T],marker=(:circle,3,.6,:red,stroke(0,.2,:black)),background_color=RGB(0.2, 0.2, 0.2),legend=false)
title!("Reduced System for Chimera States")
xaxis!(L"t",font(15))
yaxis!("Input",font(15))
