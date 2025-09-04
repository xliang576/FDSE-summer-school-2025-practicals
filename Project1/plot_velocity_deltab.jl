
Δb = [0.1,0.2,0.5,1,2,4,6] 
slope_u = [0.15,0.26,0.39,0.50,0.77,1.11,1.39]


p = plot(log.(Δb),log.(slope_u),  xlabel = "log(Δb)", ylabel = "log(slope)")
plot!(log.(Δb), log.(sqrt.(0.5*Δb)))

savefig(p,"log.png")