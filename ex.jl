using Test
import DPFEHM
import JLD
import Zygote

include("setup_model.jl")
neighbors, areasoverlengths, Ks, Qs, dirichletnodes, dirichleths, xs, ys, zs, topnodes = JLD.load("waffledata.jld", "neighbors", "areasoverlengths", "hycos", "sources", "dirichletnodes", "dirichletheads", "xs", "ys", "zs", "topnodes")
@time h = DPFEHM.groundwater_steadystate(Ks, neighbors, areasoverlengths, dirichletnodes, dirichleths, Qs; reltol=1e-12)
obsnodes = rand(1:length(Qs), 30)
h_obs = h[obsnodes] + 0.03 * randn(length(obsnodes))
function objfunc(logKs)
	h = DPFEHM.groundwater_steadystate(exp.(logKs), neighbors, areasoverlengths, dirichletnodes, dirichleths, Qs; reltol=1e-12)
	return sum((h_obs - h[obsnodes]) .^ 2)
end
@time gradient = Zygote.gradient(objfunc, log.(Ks))[1]

function checkgradientquickly(f, x0, gradf, n; delta::Float64=1e-8, kwargs...)
	indicestocheck = sort(collect(1:length(x0)), by=i->abs(gradf[i]), rev=true)[1:n]
	f0 = f(x0)
	for i in indicestocheck
		x = copy(x0)
		x[i] += delta
		fval = f(x)
		grad_f_i = (fval - f0) / delta
		@test isapprox(gradf[i], grad_f_i; kwargs...)
	end
end

checkgradientquickly(objfunc, log.(Ks), gradient, 3; delta=1e-6, rtol=1e-2)
