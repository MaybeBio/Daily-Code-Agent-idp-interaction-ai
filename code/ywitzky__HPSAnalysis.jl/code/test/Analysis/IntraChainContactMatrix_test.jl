
Sim = HPSAnalysis.SimData()

Sim.NAtoms=20
Sim.NChains=2
Sim.NSteps=5
Sim.ChainStart=[1, 13]
Sim.ChainStop=[12, 20]
Sim.ChainLength=[12, 8]
Sim.EquilibrationTime=1
Sim.RGMeasureStep=1

Sim.x_uw = zeros(Sim.NAtoms, Sim.NSteps)
Sim.y_uw = zeros(Sim.NAtoms, Sim.NSteps)
Sim.z_uw = zeros(Sim.NAtoms, Sim.NSteps)

for i in 1:5
    Sim.x_uw[1:20, i] .= 1:20
    Sim.y_uw[1:20, i] .= 1:20
    Sim.z_uw[1:20, i] .= 1:20
end

Result = HPSAnalysis.computeIntraChainContactMatrix(Sim)

res1 = zeros(Float32, 12,12)
res2 = zeros(Float32, 8,8)
for step in 1:5
    for i in 1:12
        for j in 1:12
            dist =  (Sim.x_uw[i,step]-Sim.x_uw[j,step])^2
            dist += (Sim.y_uw[i,step]-Sim.y_uw[j,step])^2
            dist += (Sim.z_uw[i,step]-Sim.z_uw[j,step])^2
            res1[i,j] += sqrt(dist)
        end
    end
end
res1 ./=5
for i in 1:12
    for j in 1:12
        if i != j
            res1[i,j] = log.(res1[i,j])
        end
    end
end


for step in 1:5
    for (i_rel,i) in enumerate(13:20)
        for (j_rel,j) in enumerate(13:20)
            dist =  (Sim.x_uw[i,step]-Sim.x_uw[j,step])^2
            dist += (Sim.y_uw[i,step]-Sim.y_uw[j,step])^2
            dist += (Sim.z_uw[i,step]-Sim.z_uw[j,step])^2
            res2[i_rel, j_rel] += sqrt(dist)
        end
    end
end
res2 ./=5
for i in 1:8
    for j in 1:8
        if i != j
            res2[i,j] = log.(res2[i,j])
        end
    end
end


@testset "Intra chain contact map" begin
    @test Result == [res1, res2]
end 