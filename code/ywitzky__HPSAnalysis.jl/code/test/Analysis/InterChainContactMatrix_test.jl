Sim = HPSAnalysis.SimData()

Sim.NAtoms=50
Sim.NChains=5
Sim.NSteps=200
Sim.ChainStart =[1 ,11,21,31,41]
Sim.ChainStop = [10,20,30,40,50]
Sim.ChainLength=[10,10,10,10,10]
Sim.EquilibrationTime=1
Sim.RGMeasureStep=1
Sim.IDs = vcat(collect(1:10), collect(1:10), collect(1:10), collect(1:10), collect(1:10))
Sim.Masses = ones(50)
Sim.ChainMasses = [10.0,10.0,10.0,10.0,10.0]

Sim.x_uw = rand(Sim.NAtoms, Sim.NSteps)*1.4.-0.7
Sim.y_uw = rand(Sim.NAtoms, Sim.NSteps)*1.2.-0.6
Sim.z_uw = rand(Sim.NAtoms, Sim.NSteps).-0.5
Sim.BoxSize = transpose([-0.7; 0.7;; -0.6; 0.6;; -0.5; 0.5]) 


Sim.x = Sim.x_uw
Sim.y = Sim.y_uw
Sim.z = Sim.z_uw

Sim.BoxLength = [1.4,1.2,1.0]


HPSAnalysis.computeCOM!(Sim)
HPSAnalysis.computeClustersByChainCOM(Sim; Cutoff=2.0) ### all proteins will be in the same cluster

### check if naiv implementation matches the quick solution
CutoffDict = Dict(i=> 0.05+i/100.0 for i in 1:10) 


### very naiv solutions
matrix = zeros(10,10)
submatrix = [zeros(10,10) for _ in 1:5]

for t in 1:Sim.NSteps
    for (I, ioff) in enumerate([0,10,20,30,40])
        for (J,joff) in enumerate([0,10,20,30,40])
            if ioff==joff continue end
            for i in 1:10
                for j in 1:10
                    x_dist =  (Sim.x[ioff+i,t] - Sim.x[joff+j,t])#%Sim.BoxLength[1])
                    y_dist =  (Sim.y[ioff+i,t] - Sim.y[joff+j,t])#%Sim.BoxLength[2])
                    z_dist =  (Sim.z[ioff+i,t] - Sim.z[joff+j,t])#%Sim.BoxLength[3])

                    x_dist -= Sim.BoxLength[1]*round(Int32, x_dist/Sim.BoxLength[1])
                    y_dist -= Sim.BoxLength[2]*round(Int32, y_dist/Sim.BoxLength[2])
                    z_dist -= Sim.BoxLength[3]*round(Int32, z_dist/Sim.BoxLength[3])

                    dist = sqrt(x_dist^2+y_dist^2+z_dist^2)
                    if dist < 1.1*(CutoffDict[i]+CutoffDict[j])/2.0
                        matrix[i,j] += 1
                        submatrix[I][i,j] += 1
                    end
                end
            end
        end
    end
end
matrix ./= 5.0*Sim.NSteps

### somewhat naiv solution for tests with realistic system sizes
HPSAnalysis.computeInterChainContactMatrix_Naiv(Sim; CutoffDict=CutoffDict, rel_cutoff=1.1)
naiv_solution = deepcopy(Sim.ContactMatrices)

HPSAnalysis.computeInterChainContactMatrix(Sim;CutoffDict=CutoffDict, rel_cutoff=1.1)
fast_solution = deepcopy(Sim.ContactMatrices)
fast_error = deepcopy(Sim.ContactMatricesError)

@test all(matrix .≈ naiv_solution)
@test all(naiv_solution .≈ fast_solution)


error = zeros(10,10)
for i in 1:5
    submatrix[i] .+= transpose(submatrix[i]) ### symmetrize matrix
    submatrix[i] ./= Sim.NSteps*2
    error .+= @.  (submatrix[i] - matrix)^2
end

error = sqrt.(error)./sqrt(5)
@test all(error .≈ fast_error)

### this version isnt really faster, but it was thought to be
HPSAnalysis.computeInterChainContactMatrix(Sim;CutoffDict=CutoffDict, rel_cutoff=1.1)
faster_solution = deepcopy(Sim.ContactMatrices)
faster_error = deepcopy(Sim.ContactMatricesError)

@test all(naiv_solution .≈ faster_solution)
@test all(error .≈ faster_error)


### check if positions which are exactly at the BoxLength and are therefore not wrapped
### get a hash which is still "inside" the box
CL_sizes = extrema.(axes(Sim.CellList))### get the first hashes outside the box
Sim.x_uw[1,1]=Sim.BoxSize[1,1]
Sim.x_uw[2,2]=Sim.BoxSize[1,2]
Sim.y_uw[3,3]=Sim.BoxSize[2,1]
Sim.y_uw[4,4]=Sim.BoxSize[2,2]
Sim.x_uw[5,5]=Sim.BoxSize[3,1]
Sim.x_uw[6,6]=Sim.BoxSize[3,2]

res = []
for step in 1:6
    Sim.CellStep[1] = step
    tuples = HPSAnalysis.getUniqueCellsOfChain(Sim, 1)
    x_ex =  extrema(getindex.(tuples, 1))
    y_ex =  extrema(getindex.(tuples, 2))
    z_ex =  extrema(getindex.(tuples, 3))

    ### true if any of hashes is "outside" the box
    val= any(map(x->x∈CL_sizes[1],x_ex)) || any(map(x->x∈CL_sizes[2],y_ex)) || any(map(x->x∈CL_sizes[3],z_ex))
    push!(res, val)
end

@test all( .!(res)) ### 
### get getUniqueCellsOfChain test with x=BoxLength



