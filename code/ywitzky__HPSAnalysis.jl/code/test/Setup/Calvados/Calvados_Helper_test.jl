using HPSAnalysis,Test
using HPSAnalysis.BioData
using GSDFormat

SimulationType_list=["HPS-Alpha","else"]
PH_list=[7.5,8.0]

@testset "DetermineAtomTypes" begin
    for pH in PH_list
        for SimulationType in SimulationType_list
            Sequences=["MRPVFV","MRPVF","MRPV","MRP"]
            AtomTypes_test=Set(join(Sequences))
            AaToId_test = Dict{Char,Int32}()
            for (index, value) in enumerate(AtomTypes_test)
                AaToId_test[value]=index
            end
            Long=Set()
            OneToCharge_test=Dict() 
            OneToMass_test=deepcopy(BioData.AaToWeight)
            OneToSigma_test=Dict() 
            OneToLambda_test=Dict() 
            OneToHPSDihedral0110_test=deepcopy(BioData.OneToHPSDihedral0110)
            OneToHPSDihedral1001_test=deepcopy(BioData.OneToHPSDihedral1001)
            LongAtomTypesToRes_test=Dict{Char,Tuple{Char,Bool}}()

            if SimulationType=="HPS-Alpha"
                OneToCharge_test = deepcopy(BioData.OneToHPSCharge)
                OneToLambda_test = deepcopy(BioData.OneToHPSUrryLambda)
                OneToSigma_test  = deepcopy(BioData.OneToHPSCalvadosSigma)
            else
                OneToCharge_test = deepcopy(BioData.OneToHPSCharge)
                OneToLambda_test = deepcopy(BioData.OneToCalvados2Lambda)
                OneToSigma_test = deepcopy(BioData.OneToHPSCalvadosSigma)
            end

            LongAtomTypes_test=union(AtomTypes_test, Long)
            IdToAa_test=Dict((v=>k) for (k,v) in AaToId_test)
            ResToLongAtomType_test = Dict{Any, Any}()

            (AtomTypes, LongAtomTypes, AaToId, IdToAa,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, OneToHPSDihedral0110, OneToHPSDihedral1001)=HPSAnalysis.Setup.DetermineAtomTypes(Sequences,SimulationType,pH)

            @test (AtomTypes_test==AtomTypes)
            @test (LongAtomTypes_test==LongAtomTypes)
            @test (AaToId_test==AaToId)
            @test (IdToAa_test==IdToAa)
            @test (ResToLongAtomType_test==ResToLongAtomType)
            @test (LongAtomTypesToRes_test==LongAtomTypesToRes)
            @test (OneToCharge_test==OneToCharge)
            @test (OneToMass_test==OneToMass)
            @test (OneToSigma_test==OneToSigma)
            @test (OneToLambda_test==OneToLambda)
            @test (OneToHPSDihedral0110_test==OneToHPSDihedral0110)
            @test (OneToHPSDihedral1001_test==OneToHPSDihedral1001)
        end
    end
end

PH_list=[7.5,8.0]
@testset "CalvadosSetup" begin
    for pH in PH_list
        Sequences=["MRPVFV","MRPVF","MRPV","MRP"]
        AtomTypes_test=Set(join(Sequences))
        AaToId_test = Dict{Char,Int32}()
        for (index, value) in enumerate(AtomTypes_test)
            AaToId_test[value]=index
        end
        Long=Set()
        OneToCharge_test=Dict() 
        OneToMass_test=deepcopy(BioData.AaToWeight)
        OneToSigma_test=Dict() 
        OneToLambda_test=Dict() 
        OneToHPSDihedral0110_test=deepcopy(BioData.OneToHPSDihedral0110)
        OneToHPSDihedral1001_test=deepcopy(BioData.OneToHPSDihedral1001)
        LongAtomTypesToRes_test=Dict{Char,Tuple{Char,Bool}}()

        AaToId_test['a']=Int32(6)
        AaToId_test['b']=Int32(7)
        AaToId_test['c']=Int32(8)
        AaToId_test['d']=Int32(9)
        Long=Set(['a','b','c','d'])
        LongAtomTypesToRes_test['a']=('M',1)
        LongAtomTypesToRes_test['b']=('V',0)
        LongAtomTypesToRes_test['c']=('F',0)
        LongAtomTypesToRes_test['d']=('P',0)

        LongAtomTypes_test=union(AtomTypes_test, Long)

        OneToCharge_test = deepcopy(BioData.OneToHPSCharge)
        OneToLambdaC2_test = deepcopy(BioData.OneToCalvados2Lambda)
        OneToLambdaC3_test = deepcopy(BioData.OneToCalvados3Lambda)
        OneToSigma_test  = deepcopy(BioData.OneToHPSCalvadosSigma)
        OneToCharge_test['H'] = 1. / ( 1 + 10^(pH-6) ) 

        OneToCharge_test['a']=OneToCharge_test['M']+1
        OneToMass_test['a']=OneToMass_test['M']+2.0
        OneToSigma_test['a']=OneToSigma_test['M']
        OneToLambdaC2_test['a']=OneToLambdaC2_test['M']
        OneToLambdaC3_test['a']=OneToLambdaC3_test['M']
        OneToHPSDihedral0110_test['a']=OneToHPSDihedral0110_test['M']
        OneToHPSDihedral1001_test['a']=OneToHPSDihedral1001_test['M']

        OneToCharge_test['b']=OneToCharge_test['V']-1
        OneToMass_test['b']=OneToMass_test['V']+16.0
        OneToSigma_test['b']=OneToSigma_test['V']
        OneToLambdaC2_test['b']=OneToLambdaC2_test['V']
        OneToLambdaC3_test['b']=OneToLambdaC3_test['V']
        OneToHPSDihedral0110_test['b']=OneToHPSDihedral0110_test['V']
        OneToHPSDihedral1001_test['b']=OneToHPSDihedral1001_test['V']

        OneToCharge_test['c']=OneToCharge_test['F']-1
        OneToMass_test['c']=OneToMass_test['F']+16.0
        OneToSigma_test['c']=OneToSigma_test['F']
        OneToLambdaC2_test['c']=OneToLambdaC2_test['F']
        OneToLambdaC3_test['c']=OneToLambdaC3_test['F']
        OneToHPSDihedral0110_test['c']=OneToHPSDihedral0110_test['F']
        OneToHPSDihedral1001_test['c']=OneToHPSDihedral1001_test['F']

        OneToCharge_test['d']=OneToCharge_test['P']-1
        OneToMass_test['d']=OneToMass_test['P']+16.0
        OneToSigma_test['d']=OneToSigma_test['P']
        OneToLambdaC2_test['d']=OneToLambdaC2_test['P']
        OneToLambdaC3_test['d']=OneToLambdaC3_test['P']
        OneToHPSDihedral0110_test['d']=OneToHPSDihedral0110_test['P']
        OneToHPSDihedral1001_test['d']=OneToHPSDihedral1001_test['P']

        IdToAa_test=Dict((v=>k) for (k,v) in AaToId_test)
        ResToLongAtomType_test = Dict{Any, Any}(('F', false) => 'c', ('M', true) => 'a', ('V', false) => 'b', ('P', false) => 'd')

        (AtomTypes_C2, LongAtomTypes_C2, AaToId_C2, IdToAa_C2,ResToLongAtomType_C2, LongAtomTypesToRes_C2, OneToCharge_C2, OneToMass_C2, OneToSigma_C2, OneToLambda_C2, OneToHPSDihedral0110_C2, OneToHPSDihedral1001_C2)=HPSAnalysis.Setup.DetermineAtomTypes(Sequences,"Calvados2",pH)
        Sequences=["MRPVFV","MRPVF","MRPV","MRP"]
        (AtomTypes_C3, LongAtomTypes_C3, AaToId_C3, IdToAa_C3,ResToLongAtomType_C3, LongAtomTypesToRes_C3, OneToCharge_C3, OneToMass_C3, OneToSigma_C3, OneToLambda_C3, OneToHPSDihedral0110_C3, OneToHPSDihedral1001_C3)=HPSAnalysis.Setup.DetermineAtomTypes(Sequences,"Calvados3",pH)

        @test (AtomTypes_test==AtomTypes_C2 && AtomTypes_test==AtomTypes_C3)
        @test (LongAtomTypes_test==LongAtomTypes_C2&&LongAtomTypes_test==LongAtomTypes_C3)
        @test (AaToId_test==AaToId_C2&&AaToId_test==AaToId_C3)
        @test (IdToAa_test==IdToAa_C2&&IdToAa_test==IdToAa_C3)
        @test (ResToLongAtomType_test==ResToLongAtomType_C2)
        @test (ResToLongAtomType_test==ResToLongAtomType_C3)
        @test (LongAtomTypesToRes_test==LongAtomTypesToRes_C2&&LongAtomTypesToRes_test==LongAtomTypesToRes_C3)
        @test (OneToCharge_test==OneToCharge_C2&&OneToCharge_test==OneToCharge_C3)
        @test (OneToMass_test==OneToMass_C2&&OneToMass_test==OneToMass_C3)
        @test (OneToSigma_test==OneToSigma_C2&&OneToSigma_test==OneToSigma_C3)
        @test (OneToLambdaC2_test==OneToLambda_C2&&OneToLambdaC3_test==OneToLambda_C3)
        @test (OneToHPSDihedral0110_test==OneToHPSDihedral0110_C2&&OneToHPSDihedral0110_test==OneToHPSDihedral0110_C3)
        @test (OneToHPSDihedral1001_test==OneToHPSDihedral1001_C2&&OneToHPSDihedral1001_test==OneToHPSDihedral1001_C3)
    end
end

@testset "Yukawa Interaction" begin
    SimulationType_list=["Calvados2","Calvados3","else"]
    SaltConcentration_list=[0.5,0.75]
    Temperature_list=[300,250]
    for Temperature in Temperature_list
        for SimulationType in SimulationType_list
            for SaltConcentration in SaltConcentration_list
                if SimulationType in ["Calvados2","Calvados3"]
                    e = 1.6021766### Charge of electron
                    e_0 = 8.854188### vacuum permitivity
                    NA = 6.022#14086# 1/mol Avogadro constant
                    kb = 0.00831446262
                    kT = kb*Temperature

                    epsilon_r=5321.0/Temperature+233.76-0.9297*Temperature+1.417*10.0^(-3)*Temperature^2-8.292*10.0^(-7)*Temperature^3
                    lamb = e^2/(4.0*pi*e_0*epsilon_r)*NA*1000/kT
                    kappa=sqrt((8.0*pi*lamb*SaltConcentration*NA/10))

                    ϵ_r_test=epsilon_r
                    κ_test=kappa
                else
                    ϵ_r_test=80.0
                    κ_test=10.0
                end

                ϵ_r, κ=HPSAnalysis.Setup.DetermineYukawaInteractions(;SimulationType,Temperature,SaltConcentration)

                @test isapprox(ϵ_r_test,ϵ_r; atol=1e-5)
                @test isapprox(κ_test,κ; atol=1e-5)
            end
        end
    end
end
