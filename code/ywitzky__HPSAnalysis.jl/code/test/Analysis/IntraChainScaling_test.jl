

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

for step in 1:Sim.NSteps
    for i in Sim.ChainStart[1]:Sim.ChainStop[1]
        Sim.x_uw[i, step] = i
        Sim.y_uw[i, step] = i
        Sim.z_uw[i, step] = i
    end

    for i in Sim.ChainStart[2]:Sim.ChainStop[2]
        Sim.x_uw[i, step] = 20.0-i+step
        Sim.y_uw[i, step] = 20.0-i+step
        Sim.z_uw[i, step] = 0
    end
end

Results = HPSAnalysis.computeIntraChainScalingNaiv(Sim)

res1= [ i<j ? (i-j)^2 : 0  for i in 1:12,  j in 1:12].*3
res2= [ i<j ? (i-j)^2 : 0  for i in 1:8,  j in 1:8].*2

@test Results == [res1, res2]

Results = HPSAnalysis.computeIntraChainScalingSlidingWindow(Sim)

res1= collect([ 3*i^2   for i in 1:11])
res2= collect([ 2*i^2  for i in 1:7])

@test Results == [res1, res2]


for step in 1:Sim.NSteps
    for i in Sim.ChainStart[1]:Sim.ChainStop[1]
        Sim.x_uw[i, step] = i
        Sim.y_uw[i, step] = i >5 ? i+5 : i
        Sim.z_uw[i, step] = i
    end

    for i in Sim.ChainStart[2]:Sim.ChainStop[2]
        Sim.x_uw[i, step] = 20.0-i+step
        Sim.y_uw[i, step] = step >=3 ? 20.0-i*3+step : 20.0-i+step
        Sim.z_uw[i, step] = 0
    end
end

Results = HPSAnalysis.computeIntraChainScalingNaiv(Sim)

res1= Float32.([ i<j ?  2*(i-j)^2+ (Sim.y_uw[i,1]-Sim.y_uw[j,1])^2 : 0  for i in 1:12,  j in 1:12])
res2= [ i<j ? (i-j)^2 : 0  for i in 1:8, j in 1:8].*7  .+  [ i<j ? (3*i-3*j)^2 : 0 for i in 1:8, j in 1:8].*3
res2 = Float32.(res2)
res2./=5.0


@test Results == [res1, res2]


Results = HPSAnalysis.computeIntraChainScalingSlidingWindow(Sim)

res1= Float32.([2*i^2*(12-i)  for i in 1:11])
for i in 1:12
    for j in i+1:12
        dt = j-i
        res1[dt] += (Sim.y_uw[i, 1]-Sim.y_uw[j, 1])^2
    end
end
for i in 1:11
    res1[i] /= 12 -i 
end


res2= [7*(i^2)  for i in 1:7] .+ [(3*i)^2  for i in 1:7].*3
res2 = Float32.(res2./5)

@test Results == [res1, res2]

