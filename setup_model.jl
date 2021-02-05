import DelimitedFiles
import JLD
import FEHM

if !isfile("waffledata.jld")
	if !isdir("waffle2017")#if it isn't already there, clone the repo with the mesh, etc
		run(`git clone git@gitlab.com:LANL-EM/waffle2017.git`)
	end
	if !isdir("otherdata")
		run(`tar xzf otherdata.tar.gz`)
	end

	function fehmhyco2fvhyco(xs, ys, zs, kxs, kys, kzs, neighbors)
		ks = Array{eltype(kxs)}(undef, length(neighbors))
		for i = 1:length(neighbors)
			node1, node2 = neighbors[i]
			dx = xs[node1] - xs[node2]
			dy = ys[node1] - ys[node2]
			dz = zs[node1] - zs[node2]
			dist = sqrt(dx^2 + dy^2 + dz^2)
			k1 = (abs(kxs[node1] * dx) + abs(kys[node1] * dy) + abs(kzs[node1] * dz)) / dist
			k2 = (abs(kxs[node2] * dx) + abs(kys[node2] * dy) + abs(kzs[node2] * dz)) / dist
			ks[i] = sqrt(k1 * k2)
		end
		return ks
	end

	westzonenum, westnodes = FEHM.parsezone(joinpath(dirname(@__FILE__), "waffle2017/smoothgrid2/out_west.zonn"))
	eastzonenum, eastnodes = FEHM.parsezone(joinpath(dirname(@__FILE__), "waffle2017/smoothgrid2/out_east.zonn"))
	topzonenum, topnodes = FEHM.parsezone(joinpath(dirname(@__FILE__), "waffle2017/smoothgrid2/out_top.zonn"))
	wellzonenums, wellzonenodes = FEHM.parsezone(joinpath(dirname(@__FILE__), "waffle2017/smoothgrid2/well_screens.zonn"))
	isanode, zoneornodenums, skds, eflows, aipeds = FEHM.parseflow(joinpath(dirname(@__FILE__), "otherdata/wl.flow"))
	hycoisanode, hycozoneornodenums, kxs, kys, kzs = FEHM.parsehyco(joinpath(dirname(@__FILE__), "otherdata/w01.hyco"))
	#xs, ys, zs, _ = FEHM.parsegeo(joinpath(dirname(@__FILE__), "waffle2017/smoothgrid2/tet.geo"), false)
	coords = DelimitedFiles.readdlm(joinpath(dirname(@__FILE__), "waffle2017/smoothgrid2/tet.xyz"))'
	xs = coords[1, :]
	ys = coords[2, :]
	zs = coords[3, :]
	r36loc = [5.010628700000E+05, 5.388061300000E+05, 1.765E+03, 1.767E+03, 50]
	r36infilnodes = Int[]
	for i = 1:length(xs)
		if zs[i] >= r36loc[3] && zs[i] <= r36loc[4]
			if sqrt((xs[i] - r36loc[1])^2 + (ys[i] - r36loc[2])^2) < r36loc[end]
				push!(r36infilnodes, i)
			end
		end
	end
	zonenums = [westzonenum; eastzonenum; wellzonenums; [336000]]
	nodesinzones = [westnodes; eastnodes; wellzonenodes; [r36infilnodes]]
	flownodes, skds, eflows, aipeds = FEHM.flattenzones(zonenums, nodesinzones, isanode, zoneornodenums, skds, eflows, aipeds)
	hyconodes, kxs, kys, kzs = FEHM.flattenzones(zonenums, nodesinzones, hycoisanode, hycozoneornodenums, kxs, kys, kzs)#hyco is in m/s
	dirichletnodes = Int[]
	dirichletheads = zeros(length(xs))
	sources = zeros(length(xs))
	for i = 1:length(flownodes)
		if skds[i] > 0
			push!(dirichletnodes, flownodes[i])
			dirichletheads[flownodes[i]] = skds[i]
		else
			#this means skds[i] is a mass source in units of kg/s
			sources[flownodes[i]] = -skds[i] * 1e-3#convert of kg/s to m^3/s
		end
	end
	volumes, areasoverlengths, neighbors = FEHM.parsestor(joinpath(dirname(@__FILE__), "waffle2017/smoothgrid2/tet.stor"))
	goodindices = filter(i->neighbors[i][1] < neighbors[i][2], 1:length(neighbors))
	areasoverlengths = areasoverlengths[goodindices]
	neighbors = neighbors[goodindices]
	hycos = fehmhyco2fvhyco(xs, ys, zs, kxs, kys, kzs, neighbors)

	JLD.save("waffledata.jld", "neighbors", neighbors, "areasoverlengths", areasoverlengths, "hycos", hycos, "sources", sources, "dirichletnodes", dirichletnodes, "dirichletheads", dirichletheads, "xs", xs, "ys", ys, "zs", zs, "topnodes", topnodes)
end
