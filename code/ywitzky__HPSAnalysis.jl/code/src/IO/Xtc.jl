### based on https://github.com/wesbarnett/MolecularDynamics.l

# James W. Barnett
# jbarnet4@tulane.edu
# Julia module for reading in xtc file with libxdrfile

#=
INSTALL LIBRARY WITH:
wget ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.1.tar.gz
tar xvzf xdrfile-1.1.1.tar.gz
cd xdrfile-1.1.1
./configure --enable-shared
make -j N
sudo make install
=#


module Xtc

using Libdl

libxdrffile= Libdl.find_library("libxdrfile.so", vcat(Base.DL_LOAD_PATH, ["/usr/local/lib/"]) )

export xtc_init, read_xtc, close_xtc

mutable struct xtcType
    natoms::Int32
    step::Vector{Int32}
    time::Vector{Float32}
    box::Matrix{Float32}
    xyz::Matrix{Float32}
    prec::Array{Float32}
    xd::Ptr{Cvoid}
end

function xtc_init(xtcfile) 
    # Check if file exists
    if (~isfile(xtcfile)) 
        error(string(xtcfile," xtc file does not exist."))
    end

    # Get number of atoms in system
    natoms = Cint[0]
    stat = ccall( (:read_xtc_natoms,libxdrffile), Int32, (Ptr{UInt8}, Ptr{Cint}), xtcfile, natoms)

    # Check if we actually did open the file
    if (stat != 0)
        error(string("Failure in opening ", xtcfile))
    end

    # Get C xdrfile pointer
    xd = ccall( (:xdrfile_open,libxdrffile), Ptr{Cvoid},
        (Ptr{UInt8},Ptr{UInt8}), xtcfile,"r")

    # Assign everything to this type
    xtc = xtcType(
        natoms[],
        Cint[0],
        Cfloat[0],
        zeros(Cfloat,(3,3)),
        zeros(Cfloat,(3,convert(Int64,(natoms[])))),
        Cfloat[0],
        xd)

    return (stat, xtc)

end

function read_xtc(xtc)
    stat = ccall( (:read_xtc,libxdrffile), Int32, ( Ptr{Cvoid}, Int32,
        Ptr{Cint}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat} ), xtc.xd,
        xtc.natoms, xtc.step, xtc.time, xtc.box, xtc.xyz, xtc.prec) 

    if (stat != 0 | stat != 11)
        error("Failure in reading xtc frame.")
    end

	return stat, xtc
end

function close_xtc(xtc)

    stat = ccall( (:xdrfile_close,libxdrffile), Int32, ( Ptr{Cvoid}, ), xtc.xd)

	return stat

end

end
