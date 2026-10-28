filename_test="$(SetupTestPath)/GSD_write_test.gsd"
filename="$(SetupTestPath)/GSD_write.gsd"

using GSDFormat

@testset "writeGSDStartFile" begin
    #Not testing UseAngles=true, because it is an artefact of a previous version.
    UseAngles=false

    sequences=["MNAM"]
    N=4
    set=Set(join(sequences))
    AAToID=Dict{Char,Int}()
    max_seq_length = 4
    BoxSize=Vector{Float64}([20.0,30.0,20.0,0.0,0.0,0.0])
    coor_notreshaped=[5.0 5.0 5.0 5.0;;; 62.4 66.2 70.0 73.8;;; 0.0 0.0 0.0 0.0]
    image_notreshaped=[0 0 0 0;;; 0 0 0 0;;; 0 0 0 0]
    coor =[0.5 6.24 0.0; 0.5 6.62 0.0; 0.5 7.0 0.0; 0.5 7.38 0.0]
    image=[0 0 0; 0 0 0; 0 0 0; 0 0 0]

    for (num,atom) in enumerate(set)
        AAToID[atom]=num
    end

    snapshot=GSDFormat.Frame()
    snapshot.configuration.step=1
    snapshot.configuration.dimensions=3
    snapshot.configuration.box=BoxSize./10.0

    snapshot.particles.N=N
    snapshot.particles.position=coor
    IDToAA=Dict(value=>key for (key,value) in AAToID)
    snapshot.particles.types = ["M", "A", "N"]
    snapshot.particles.typeid = [0, 2, 1, 0]
    snapshot.particles.image = image
    mass_charge=Vector{Float64}([1.0,2.0,3.0,1.0])
    snapshot.particles.mass = mass_charge
    snapshot.particles.charge = mass_charge
    snapshot.particles.diameter = [0.1, 0.3, 0.2, 0.1]
    snapshot.particles.body = [-2,-2,-1,-1]

    snapshot.bonds.N=N-1
    snapshot.bonds.types = ["O-O"]
    snapshot.bonds.typeid = zeros(Int32, N-1)
    snapshot.bonds.group = [0 1; 1 2; 2 3]
    
    DiMap=Dict()
    DiList=Matrix{Int}(undef,4,4)

    file = GSDFormat.open(filename_test, 'w')
    GSDFormat.append(file, snapshot)
    GSDFormat.close(file)

    HPSAnalysis.Setup.writeGSDStartFile(filename,N,N-1,N-2,N-3,BoxSize,coor_notreshaped,AAToID,sequences,image_notreshaped,mass_charge,mass_charge,DiMap,DiList,AAToID,false, "test", [],[-2,-2,-1,-1], )

    @test files_are_equal(filename_test,filename)
end
rm(filename_test)
rm(filename)