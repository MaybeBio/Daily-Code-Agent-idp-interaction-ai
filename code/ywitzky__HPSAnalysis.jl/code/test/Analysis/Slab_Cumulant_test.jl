Sim = HPSAnalysis.SimData()

#Sim.NAtoms=
Sim.NChains = 1
Sim.ChainMasses =ones(Sim.NChains)
TotalMass= sum(Sim.ChainMasses)
#Sim.Charges = ones(Sim.NAtoms)
Sim.SlabAxis=2

Sim.NSteps = 1

### create one particle per unit volume in dilute range
diluterange = vcat(collect(-9.5:1:-5.5), collect(5.5:1:9.5))
denserange = -4.5:1:4.5
shortrange = -5.5:1:5.5


Sim.Clusters = [[[Int32(1)]]]
Sim.COM = zeros(Float32, 1,3,1)
Sim.ClusterRange=1:1

conv_fac=1.66053906660 ### conversion from u/AA^3 to kg/L which is done in the function



@testset "Binder Cumulant in Slab Simulations" begin
    ### Check getUnwrappedSlabCoordinate
    Sim.SlabAxis = 1
    @test Sim.x === HPSAnalysis.getSlabCoordinate(Sim; Unwrapped=false)
    @test Sim.x_uw === HPSAnalysis.getSlabCoordinate(Sim; Unwrapped=true) 

    Sim.SlabAxis = 2
    @test Sim.y === HPSAnalysis.getSlabCoordinate(Sim; Unwrapped=false)
    @test Sim.y_uw === HPSAnalysis.getSlabCoordinate(Sim; Unwrapped=true)

    Sim.SlabAxis = 3
    @test Sim.z === HPSAnalysis.getSlabCoordinate(Sim; Unwrapped=false)
    @test Sim.z_uw === HPSAnalysis.getSlabCoordinate(Sim; Unwrapped=true)

    ### Check getVoxelIndex
    @test  3 == HPSAnalysis.getVoxelIndex(11.0, 5.0, 0)
    @test -2 == HPSAnalysis.getVoxelIndex(-10.0, 5.0, 0)
    @test  0 == HPSAnalysis.getVoxelIndex(-10.0, 5.0, 2)
    @test -1 == HPSAnalysis.getVoxelIndex(-10.0, 3.0, 2)

    ### check centering according to Cluster COM
    Test_pos = [1.0; 2.0; 3.0; 4.0; 5.0;; 6.0; 7.0; 8.0; 9.0; 10.0;]
    Len=100.0

    @test all( (Test_pos.- 5.0)[:] .== [HPSAnalysis.getRecenteredPositions(Test_pos, atom,1 , step, [5.0], Len, 1.0/Len) for step in 1:2 for atom in 1:5])

    Len=5.0
    @test all( [1.0,2.0,-2.0,-1.0,0.0,1.0,2.0,-2.0,-1.0,0.0] .≈ [HPSAnalysis.getRecenteredPositions(Test_pos, atom,1 , step, [0.0], Len, 1.0/Len) for step in 1:2 for atom in 1:5  ])

    @test all( [-1.0,0.0,1.0,2.0,-2.0,-1.0,0.0,1.0,2.0,-2.0] .≈ [HPSAnalysis.getRecenteredPositions(Test_pos, atom,1 , step, [2.0], Len, 1.0/Len) for step in 1:2 for atom in 1:5])

    ### testing the core function

    ### create a slab in each possible Slab Axis
    for Sim.SlabAxis in 1:3
        ### start with constant densities
        xrange_dilute= Sim.SlabAxis==1 ? diluterange : shortrange
        xrange_dense = Sim.SlabAxis==1 ? denserange : shortrange 

        yrange_dilute= Sim.SlabAxis==2 ? diluterange : shortrange
        yrange_dense = Sim.SlabAxis==2 ? denserange : shortrange 

        zrange_dilute= Sim.SlabAxis==3 ? diluterange : shortrange
        zrange_dense = Sim.SlabAxis==3 ? denserange : shortrange 
        
        x_dilute = [x for x in xrange_dilute for y in yrange_dilute for z in zrange_dilute ]
        y_dilute = [y for x in xrange_dilute for y in yrange_dilute for z in zrange_dilute ]
        z_dilute = [z for x in xrange_dilute for y in yrange_dilute for z in zrange_dilute ]

        x_dense_const = [x for _ in 1:5 for x in xrange_dense for y in yrange_dense for z in zrange_dense ]
        y_dense_const = [y for _ in 1:5 for x in xrange_dense for y in yrange_dense for z in zrange_dense ]
        z_dense_const = [z for _ in 1:5 for x in xrange_dense for y in yrange_dense for z in zrange_dense ]

        Sim.BoxSize= [-6.0 6.0; -6.0 6.0; -6.0 6.0]
        Sim.BoxSize[Sim.SlabAxis,1]=-10.0
        Sim.BoxSize[Sim.SlabAxis,2]=10.0

        Sim.BoxLength=[12.0,12.0,12.0]
        Sim.BoxLength[Sim.SlabAxis]=20.0


        Ndilute = length(x_dilute)
        Ndense = length(x_dense_const)

        Sim.NAtoms=Ndilute+Ndense
        Sim.x = zeros(Float32, Sim.NAtoms, Sim.NSteps)
        Sim.y = zeros(Float32, Sim.NAtoms, Sim.NSteps)
        Sim.z = zeros(Float32, Sim.NAtoms, Sim.NSteps)

        Sim.x[1:Ndilute,1] = x_dilute
        Sim.y[1:Ndilute,1] = y_dilute
        Sim.z[1:Ndilute,1] = z_dilute
        Sim.x[Ndilute+1:end,1] = x_dense_const
        Sim.y[Ndilute+1:end,1] = y_dense_const
        Sim.z[Ndilute+1:end,1] = z_dense_const

        Sim.x_uw = deepcopy(Sim.x)
        Sim.y_uw = deepcopy(Sim.y)
        Sim.z_uw = deepcopy(Sim.z)
        Sim.Masses = ones(Float32, Sim.NAtoms)

        dense, dilute = HPSAnalysis.computeBinderCumulantsSubBoxes(Sim, 5, 5)

        @test all(conv_fac*1.0*ones(Float32, 12, 12,1) .≈ dilute)
        @test all(conv_fac*5.0*ones(Float32, 12, 12,1) .≈ dense)

        dense, dilute = HPSAnalysis.computeBinderCumulantsSubBoxes(Sim, 4, 4)

        @test all(isapprox.(conv_fac*20.0/12.0*ones(Float32, 12, 12,1), dilute; atol=10^-6))
        @test all(conv_fac*5.0*ones(Float32, 12, 12,1) .≈ dense)

        dense, dilute = HPSAnalysis.computeBinderCumulantsSubBoxes(Sim, 6, 6)

        @test all(conv_fac*1.0*ones(Float32, 12, 12,1) .≈ dilute)
        @test all(isapprox.(conv_fac*52.0/12.0*ones(Float32, 12, 12,1), dense; atol=10^-6))
        #@test all(conv_fac*52.0/12.0*ones(Float32, 12, 12,1) .≈ dense)


        ### now impose a slope along the slab axis direction =>  density should be the same for the different cumulant voxels
        slope_(SlabAxis,x,y,z) = SlabAxis==1 ? x+6 : ( SlabAxis==2 ? y+6 : z+6)
        slope(SlabAxis,x,y,z)  = floor(Int32, slope_(SlabAxis,x,y,z) )

        x_dense_slope = [x for x in xrange_dense for y in yrange_dense for z in zrange_dense for _ in 1:slope(Sim.SlabAxis,x,y,z)]
        y_dense_slope = [y for x in xrange_dense for y in yrange_dense for z in zrange_dense for _ in 1:slope(Sim.SlabAxis,x,y,z)]
        z_dense_slope = [z for x in xrange_dense for y in yrange_dense for z in zrange_dense for _ in 1:slope(Sim.SlabAxis,x,y,z)]

        NSlope = length(x_dense_slope)
        Sim.NAtoms=Ndilute+Ndense+NSlope
        Sim.x = zeros(Float32, Sim.NAtoms, Sim.NSteps)
        Sim.y = zeros(Float32, Sim.NAtoms, Sim.NSteps)
        Sim.z = zeros(Float32, Sim.NAtoms, Sim.NSteps)


        Sim.x[1:Ndilute,1] = x_dilute
        Sim.y[1:Ndilute,1] = y_dilute
        Sim.z[1:Ndilute,1] = z_dilute
        Sim.x[Ndilute+1:Ndilute+Ndense,1] = x_dense_const
        Sim.y[Ndilute+1:Ndilute+Ndense,1] = y_dense_const
        Sim.z[Ndilute+1:Ndilute+Ndense,1] = z_dense_const
        Sim.x[Ndilute+Ndense+1:end,1] = x_dense_slope
        Sim.y[Ndilute+Ndense+1:end,1] = y_dense_slope
        Sim.z[Ndilute+Ndense+1:end,1] = z_dense_slope

        Sim.x_uw = deepcopy(Sim.x)
        Sim.y_uw = deepcopy(Sim.y)
        Sim.z_uw = deepcopy(Sim.z)
        Sim.Masses = ones(Float32, Sim.NAtoms)

        dense, dilute = HPSAnalysis.computeBinderCumulantsSubBoxes(Sim, 5, 5)

        @test all(conv_fac*1.0*ones(Float32, 12, 12,1) .≈ dilute)
        @test all(conv_fac*10.5*ones(Float32, 12, 12,1) .≈ dense)


        ### now impose a slope both non slab axis direction =>  density should be gradient for the different cumulant voxels
        AllAxis = [1,2,3]
        deleteat!(AllAxis, Sim.SlabAxis)

        x_dense_slope = [x for x in xrange_dense for y in yrange_dense for z in zrange_dense for _ in 1:(slope(AllAxis[1],x,y,z)+slope(AllAxis[2],x,y,z))]
        y_dense_slope = [y for x in xrange_dense for y in yrange_dense for z in zrange_dense for _ in 1:(slope(AllAxis[1],x,y,z)+slope(AllAxis[2],x,y,z))]
        z_dense_slope = [z for x in xrange_dense for y in yrange_dense for z in zrange_dense for _ in 1:(slope(AllAxis[1],x,y,z)+slope(AllAxis[2],x,y,z))]

        NSlope = length(x_dense_slope)
        Sim.NAtoms=Ndilute+Ndense+NSlope
        Sim.x = zeros(Float32, Sim.NAtoms, Sim.NSteps)
        Sim.y = zeros(Float32, Sim.NAtoms, Sim.NSteps)
        Sim.z = zeros(Float32, Sim.NAtoms, Sim.NSteps)


        Sim.x[1:Ndilute,1] = x_dilute
        Sim.y[1:Ndilute,1] = y_dilute
        Sim.z[1:Ndilute,1] = z_dilute
        Sim.x[Ndilute+1:Ndilute+Ndense,1] = x_dense_const
        Sim.y[Ndilute+1:Ndilute+Ndense,1] = y_dense_const
        Sim.z[Ndilute+1:Ndilute+Ndense,1] = z_dense_const
        Sim.x[Ndilute+Ndense+1:end,1] = x_dense_slope
        Sim.y[Ndilute+Ndense+1:end,1] = y_dense_slope
        Sim.z[Ndilute+Ndense+1:end,1] = z_dense_slope

        Sim.x_uw = deepcopy(Sim.x)
        Sim.y_uw = deepcopy(Sim.y)
        Sim.z_uw = deepcopy(Sim.z)
        Sim.Masses = ones(Float32, Sim.NAtoms)

        dense, dilute = HPSAnalysis.computeBinderCumulantsSubBoxes(Sim, 5, 5)

        const_dilute = conv_fac*1.0*ones(Float32, 12, 12,1)
        const_dense = conv_fac*5.0*ones(Float32, 12, 12,1)
        dual_slope = conv_fac*1.0*[u+v for u in 0:11, v in 0:11]

        @test all(const_dilute .≈ dilute)
        @test all(const_dense+dual_slope .≈ dense)
    end

end



