const ChoosenFloatType=Float32
const ChoosenIntType=Int32
using Revise, OffsetArrays
using StaticArrays


mutable struct SimData{FloatType<:AbstractFloat,IntType<:Integer}
    HOOMD::Bool
    SimulationName::String
    EquilibrationTime::IntType
    NSteps::IntType
    NRestarts::IntType
    NAtoms::IntType
    NBonds::IntType
    NAngles::IntType
    NDihedrals::IntType
    NChains::IntType
    NAtomTypes::IntType
    NBondTypes::IntType
    NAngleTypes::IntType
    NDihedralTypes::IntType
    TargetTemp::FloatType
    SlabAxis::IntType
    BoxSize::Array{FloatType}
    BoxLength::Array{FloatType}
    Sequences::Vector{String}
    
    ### variables mostly used for HREMD
    LammpsVariables::Dict{String,Any}
    NReplica::IntType
    EnergyDict::Dict{String, IntType}
    Energies::Matrix{FloatType}
    EnergyWriteOutFreq::IntType
    TrajWriteOutFreq::IntType
    
    xFilePath::String
    yFilePath::String
    zFilePath::String
    xio::Union{IOStream, Nothing}
    yio::Union{IOStream, Nothing}
    zio::Union{IOStream, Nothing}
    x::Matrix{FloatType} ### Atom, Step
    y::Matrix{FloatType} ### Atom, Step
    z::Matrix{FloatType} ### Atom, Step
    #xyz::Array{FloatType}
    StepFrequency::IntType
    Reduce::Rational

    x_uw_FilePath::String
    y_uw_FilePath::String
    z_uw_FilePath::String
    x_uw_io::Union{IOStream, Nothing}
    y_uw_io::Union{IOStream, Nothing}
    z_uw_io::Union{IOStream, Nothing}
    x_uw::Matrix{FloatType} ### Atom, Step
    y_uw::Matrix{FloatType} ### Atom, Step
    z_uw::Matrix{FloatType} ### Atom, Step

    BasePath::String
    PlotPath::String
    DataPath::String
    TrajectoryFile::String
    EnergyFile::String
    StartFile::String ### contains bonds, angle and dihedral definitions of lammps
    BigDataList::Array{Symbol}

    IDs::Vector{IntType}
    Charges::Vector{FloatType}
    Masses::Vector{FloatType}
    ChainMasses::Vector{FloatType}
    COP::Array{FloatType}
    COP_uw::Array{FloatType}
    COM::Array{FloatType}
    COM_uw::Array{FloatType}
    IDToResName::Dict{IntType, String}
    IDToMasses::Vector{FloatType}
    IDToSigmas::Vector{FloatType}
    IDToLambdas::Vector{FloatType}
    ChainStart::Vector{IntType}
    ChainStop::Vector{IntType}
    ChainLength::Vector{IntType}
    MaxChainLength::IntType

    RGMeasureStep::IntType
    RGMean::FloatType
    RGErr::FloatType
    RGComponentSeries::Array{FloatType}
    RGSeries::Array{FloatType}
    RGAutocorr::Array{FloatType}
    RGMeanByChain::Vector{FloatType}
    RGErrByChain::Vector{FloatType}
    REESeries::Array{FloatType}
    REEVecSeries::Array{FloatType}
    REEComponents::Array{FloatType}
    InertiaTensorEigVals::Array{FloatType}
    ShapeAsymmetry::Array{FloatType}
    ParallelInertiaTensor::Array{FloatType}
    AspectRatio::Array{FloatType}
    Asphericity::Array{FloatType}
    Acylindricity::Array{FloatType}
    MeanShapeAsymmetry::Array{FloatType}
    MeanAspectRatio::Array{FloatType}
    MeanAsphericity::Array{FloatType}
    MeanAcylindricity::Array{FloatType}

    Resolution::FloatType
    SlabHistogramSeries::OffsetArray{FloatType}
    NegativeSlabHistogramSeries::OffsetArray{FloatType}
    PositiveSlabHistogramSeries::OffsetArray{FloatType}
    WeightMapXY::OffsetArray{FloatType}
    WeightMapXZ::OffsetArray{FloatType}
    WeightMapYZ::OffsetArray{FloatType}
    OffsetInd::Array{IntType}

    BondAngles::Array{FloatType}
    AvgBondAngles::Array{FloatType}
    CosBondAngles::Array{FloatType}
    AvgCosBondAngles::Array{FloatType}
    BondAngleHist::Array{FloatType}
    LocalPersistenceLength::Array{FloatType}
    VarLocalPersistenceLength::Array{FloatType}
    TorsionAngles::Array{FloatType}
    TorsionHist::OffsetArray{FloatType}
    AlphaHelixProb::Vector{FloatType}
    CellStep::Array{IntType}
    CenterBox::Array{FloatType}
    NeighBox::Array{FloatType}
    CL_Dist::Array{FloatType}

    MaxParticlesPerCell::IntType
    CellDimensions::Array{IntType}
    CellResolution::Array{FloatType}
    CellListCounters::OffsetArray{IntType}
    CellList::OffsetArray{Vector{IntType}}
    ChainNumCellList::OffsetArray{Vector{IntType}}
    PositiveCellList::OffsetArray{Vector{IntType}}
    NegativeCellList::OffsetArray{Vector{IntType}}

    ChargeAnalysisCutoff::FloatType
    ChargeContactTimeCellList::OffsetArray{IntType}
    #ChargeContactCellList::Dict{IntType, IntType}
    ChargeContactValid::Dict{Tuple{IntType, IntType}, Bool}
    ChargeContactTime::Dict{Tuple{IntType, IntType}, IntType}
    ChargeContactTimeHist::Array{FloatType}
    PosChargeToHistID::Dict{FloatType, IntType}
    NegChargeToHistID::Dict{FloatType, IntType}
    ChargeResidueContactMatrix::Array{FloatType}
    ChargeChainContactMatrix::Array{FloatType}

    MSD::Array{FloatType}
    BondHist::Vector{FloatType}
    box_id::Array{FloatType}

    ContactMatrices::Array{FloatType}
    ContactMatricesError::Array{FloatType}
    ContactMatrixCutoffs::Array{FloatType}
    ContactMatrixTypeCutoffs::Array{FloatType}

    HydropathyDecoration::Array{FloatType}
    MeanHydropathy::Array{FloatType}

    REEHistsAxes::Array{UnitRange{Int32}}
    REEHists::Array{FloatType}
    FrameWeights::Array{FloatType}

    ClusterRange::StepRange{IntType,IntType}
    Clusters::Vector{Vector{Vector{IntType}}}

    GofRRange::Vector{FloatType}
    GofR_Resolution::FloatType
    GofR_MaxRange::FloatType
    GofR::Vector{FloatType}
    DensityHist::Vector{FloatType}

    IntraChainScalingNaiv::Vector{Matrix{FloatType}}
    IntraChainScalingSlidingWindow::Vector{Vector{FloatType}}
    IntraChainContactMatrix::Vector{Matrix{FloatType}}

    BinderSlabCumulants::Vector{FloatType}
    SlabDenseCumulantBoxes::Array{FloatType}
    SlabDiluteCumulantBoxes::Array{FloatType}

    SurfaceTension::FloatType
    SurfaceTensionError::FloatType

    MeanYukawaEnergy::Matrix{FloatType}
    MeanYukawaEnergy_perAA::Matrix{FloatType}
    MeanAshbaughHatchEnergy::Matrix{FloatType}
    MeanAshbaughHatchEnergy_perAA::Matrix{FloatType}
    MeanIntraYukawaEnergy::Matrix{FloatType}
    MeanIntraYukawaEnergy_perAA::Matrix{FloatType}
    MeanIntraAshbaughHatchEnergy::Matrix{FloatType}
    MeanIntraAshbaughHatchEnergy_perAA::Matrix{FloatType}


    function SimData()
        return new{ChoosenFloatType, ChoosenIntType}(false,"",1,0, 0,0,0,0,0,0,0,0,0,0,0.0,2,zeros(3,2), zeros(3),[""],
        
        Dict{String,Any}(),0, Dict{String,ChoosenIntType}(),zeros((0,0)),0,0,
        
        "","","",nothing,nothing,nothing,zeros((0,0)), zeros((0,0)),zeros((0,0)),1,1, 

        "","","",nothing,nothing,nothing,zeros((0,0)), zeros((0,0)),zeros((0,0)),

        "","","","","","",[],

        zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),Dict(), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), 0,
        
        0,0,0,zeros(0),zeros(0),zeros(0), zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0),zeros(0), 
        1,zeros(0),zeros(0),zeros(0), zeros(0),zeros(0),
        
        zeros(0),zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0) ,zeros(0) ,
       
        [1],zeros(0),zeros(0),zeros(0),zeros(0), zeros(0),0,zeros(3,2), zeros(3),OffsetArray(zeros(0)), OffsetArray(zeros(0)), OffsetArray(zeros(0)), OffsetArray(zeros(0)),OffsetArray(zeros(0)),
       
        1., OffsetArray(zeros(0)), Dict(), Dict(), zeros(0), Dict(), Dict(), zeros(0), zeros(0),
       
        zeros(0), zeros(0), zeros(0),zeros(0),

        zeros(0),zeros(0),zeros(0),
        
        zeros(0), zeros(0),
        
        zeros(0),zeros(0), zeros(0),

        1:1:1,  Vector{Vector{Vector{ChoosenIntType}}}(),
        
        zeros(0), 0,0,zeros(0),zeros(0),
        
        [], [],[],
        
        zeros(0), zeros(0), zeros(0),
        
        0,0,
        
        zeros((0,0)),zeros((0,0)),zeros((0,0)),zeros((0,0)),zeros((0,0)),zeros((0,0)),zeros((0,0)),zeros((0,0))) 
    end
end

