### Rewrite AlphaFold Cif InputFiles to PDB in Folder
function RewriteCifToPDB(BasePath, ProteinToCif,Proteins)
    subpath = "$(BasePath)/InitFiles/PDBFiles"
    cifpath = "$(BasePath)/InitFiles/CifFiles"
    mkpath(subpath)
    mkpath(cifpath)
    for Prot in Set(Proteins)
        CifPath = ProteinToCif[Prot]
        try
            readlines(CifPath)
        catch
            @error "There is no AlphaFold data for protein $(Prot) at $(CifPath)."
        end

        cp(CifPath, "$(cifpath)/$(Prot).cif"; force=true)
        pdb_file = open("$(subpath)/$(Prot).pdb", "w")
        for line in readlines(CifPath)
            fields = strip.(split(line))
            if !isempty(fields) && fields[1] == "ATOM"
                ID = parse(Int, fields[2])
                symbole = fields[3]
                label_atom = fields[4]
                comp = fields[6]
                asym = fields[7]
                seq  = parse(Int, fields[9])
                x = parse(Float64, fields[11])
                y = parse(Float64, fields[12])
                z = parse(Float64, fields[13])
                iso_or_equiv= parse(Float64, fields[15])

                println(pdb_file, @sprintf("ATOM  %5d %4s %-3s %s%4d    %8.3f%8.3f%8.3f  1.00  %1.2f %10s", ID, label_atom, comp, asym, seq, x, y, z, iso_or_equiv, symbole))
            end
        end
        close(pdb_file)
    end
end
