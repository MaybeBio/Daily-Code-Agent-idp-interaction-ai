using ..BioData

using BioAlignments, BioSequences, BioSymbols, DataFrames

### taken from https://pubs.aip.org/aip/jcp/article/148/12/123305/1023768/Sequence-charge-decoration-dictates-coil-globule
function SequenceChargeDecorator(charges::Vector{R}) where {R<:Real} 
    N = length(charges)
    SCD = zero(R)
    sqrt_pre = sqrt.(collect(1:N))
    for m in 2:N
        for n in 1:m-1
            SCD += charges[m]*charges[n]*sqrt_pre[m-n]
        end
    end
    return SCD/N
end

function SequenceChargeDecorator_Var1(charges::Vector{R}) where {R<:Real} 
    N = length(charges)
    mean_charge = sum(charges)/N

    SCD = zero(R)
    for m in 2:N
        if charges[m] == 0 continue end
        for n in 1:m-1
            if charges[n] == 0 continue end
            SCD += (charges[m]-mean_charge)*(charges[n]-mean_charge)*sqrt(m-n)
        end
    end
    return SCD/N
end

function SequenceChargeDecorator_Var2(charges::Vector{R}) where {R<:Real} 
    N = length(charges)
    mean_charge_sqr = (sum(charges)/N)^2
    SCD = zero(R)
    for m in 2:N
        if charges[m] == 0 continue end
        for n in 1:m-1
            if charges[n] == 0 continue end
            SCD += (charges[m]*charges[n]-mean_charge_sqr)*sqrt(m-n)
        end
    end
    return SCD/N
end


### assumes Calvados2 and pH 7 , H needs to be modified for other pH values
function ConvertSequenceToCharges(Sequence::String, Type=Float32, OneToCharge=BioData.OneToCalvados2Charge; pH=7.0)
    OneToCharge['H']=1. / ( 1 + 10^(pH-6) )

    charges=zeros(Type, length(Sequence))
    for (i,residue) in enumerate(Sequence)
        charges[i]=OneToCharge[residue]
    end
    return charges
end

### DEPRECATED!!!! need monte carlo for this
### assumes charges +- 1
function createMinMaxSequences(charges::Vector{R}; InsertEmpty=true) where {R<:Real}
    q_net = sum(charges)
    q_abs = sum(abs.(charges))
    N1 = Int32(q_abs/2 + q_net/2)
    N2 = Int32(q_abs/2 - q_net/2)

    if InsertEmpty
        max_val = zeros(length(charges))
    else
        max_val = zeros(N1+N2)
    end
    prev = q_net>1 ? -1.0 : 1.0
    max_val[1:N1] .= 1
    max_val[end-N2+1:end] .= -1.0

    q_net2 = sum(max_val)
    q_abs2 = sum(abs.(max_val))

    min_val = ones(Int32(q_abs))
    min_val[1:2:end] *= -1.0
    if q_net!=0
        if q_net>0
            min_val[1:ceil(Int32, abs(q_net)/2)] .= 1*sign(q_net)
            min_val[end-floor(Int32, abs(q_net)/2):end] .= 1*sign(q_net)
        else
            min_val[1:ceil(Int32, abs(q_net)/2)+1] .= 1*sign(q_net)
            min_val[end-floor(Int32, abs(q_net)/2)+1:end] .= 1*sign(q_net)
        end
    end
    return min_val, max_val
end


### assumes charges +- : DEPRECATED!!!
function normSCDValues(data::Vector{R}, charge::Vector{R2}; func=SequenceChargeDecorator) where {R<:Real, R2<:Real}
    #min_val, max_val = createMinMaxSequences(charge)

    q_net = sum(min_val)
    q_abs = sum(abs.(min_val))

    ### naming gets weird since they invert scale
    tmp = data .* -1.0
    minVal = min(minimum(tmp),-1.0*func(min_val))
    maxVal = max(maximum(tmp),-1.0*func(max_val))

    #("tmp $tmp, minVal: $minVal, func(min_val) $(-1.0*func(min_val)) maxVal: $maxVal")
    return abs.((tmp .-minVal)./(abs(maxVal)-minVal))
end

function AlignSequences(Seq1::String, Seq2::String)
    PhosRevDict = Dict{Char, Char}('#'=>'S', '&'=>'T', '@'=>'Y')
    phosSites = [i for (i, res) in enumerate(Seq2) if res in ['#','&','@']]
    Seq2 = [ (res in ['#','&','@'] ? PhosRevDict[res] : res) for res in Seq2 ]

    Seq1_AA = convert.(AminoAcid, Vector{Char}(Seq1))
    Seq2_AA = convert.(AminoAcid, Vector{Char}(Seq2))

    scoremodel = AffineGapScoreModel(BLOSUM62, gap_open=-5, gap_extend=-1);
    result = pairalign(GlobalAlignment(), Seq1_AA, Seq2_AA, scoremodel)
    seq, ref, mat = print_pairwise_alignment_YW(alignment(result), phosSites)

end

### https://github.com/BioJulia/BioAlignments.jl/blob/ceb3fe6ce85c794adde877b0bd048905557dc1e2/src/pairwise/alignment.jl
# modified from BioAlignments/src/pairwise/alignment.jl
function print_pairwise_alignment_YW(aln::PairwiseAlignment, phosSites::Vector{I}; width::Integer=60) where{ I<:Integer}
    PhosDict = Dict{Char, Char}('S'=>'#', 'T'=>'&', 'Y'=>'@')
    seq = aln.a.seq
    ref = aln.b
    anchors = aln.a.aln.anchors
    # width of position numbers
    posw = ndigits(max(anchors[end].seqpos, anchors[end].refpos)) + 1

    i = 0
    seqpos = anchors[1].seqpos
    refpos = anchors[1].refpos
    seq_cnt = 0
    ref_cnt = 0
    seq=""
    ref=""
    mat=""
    next_xy = iterate(aln)
    while next_xy !== nothing
        (x, y), s = next_xy
        next_xy = iterate(aln ,s)

        seq_cnt += y != AA_Gap
        ref_cnt += x != AA_Gap
        seq *= seq_cnt in phosSites ? PhosDict[Char(x)] :  Char(y)

        ref *= Char(x)
        mat *= x == y ? '|' : ' '
    end

    return seq,  ref, mat
end


function ColorSequencesHelp(Sequence::String, OneToCharge=BioData.OneToCalvados2Charge)
    cur_string = ""
    prev_sign =0
    texts = Vector{String}()
    signs = Vector{Int32}()
    for res in Sequence
        q = get(OneToCharge, res,0)
        if sign(prev_sign)==sign(q)
            cur_string *= res
        else
            push!(texts, cur_string)
            push!(signs, prev_sign)
            cur_string="$(res)"
            prev_sign=sign(q)
        end
    end
    push!(texts, cur_string)
    push!(signs, prev_sign)
    return texts, signs
end

### needs to be used with python package openpyxl using PyCall
# tmp = createTextBlockString.(ColorSequencesForExcel(Sequences))
# py""" wb = openpyxl.Workbook() 
#   ws = wb.active
#   ws["A1"] = $$tmp
#   wb.save(filename)"""
function createTextBlockString(signs, texts)
    ColorDict3 = Dict(1=>4, -1=>2 , 0=>0)

    out = ""
    for (sign, text) in zip(signs, texts)
        out *= ", openpyxl.cell.rich_text.TextBlock(openpyxl.cell.text.InlineFont(rFont=\'Cambria\', color=openpyxl.styles.colors.Color(indexed=$(ColorDict3[sign]))), text=\"$text\")"
    end
    return "openpyxl.cell.rich_text.CellRichText($(out[2:end]))"
end

function ColorSequencesForExcel(Sequence::String, OneToCharge=BioData.OneToCalvados2Charge)
    texts, signs = ColorSequencesHelp(Sequence, OneToCharge)
    return createTextBlockString(signs, texts)
end


function computeKappa_Helper(charges::Vector{R}, charge_asym; width=5) where {R<: Real}
    ### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3746876/pdf/pnas.201304749.pdf
    ### Conformations of intrinsically disordered proteins are influenced by linear sequence istributions of oppositely charged residues

    N = length(charges)
    if N<width
        return -1000
    end
    sigma = zeros(R, N-width+1)
    for i in 1:N-width+1
        q_pos =  zero(R)
        N_pos = zero(Int32)
        q_neg =  zero(R)
        N_neg = zero(Int32)  
        for j in i:i+width-1
            q = charges[j]
            if q >0 
                q_pos += q 
                N_pos += 1
            elseif q<0 
                q_neg += q 
                N_neg += 1
            end
        end
        f_plus =  N_pos/width
        f_neg  =  N_neg/width

        if f_plus+f_neg ==0
            sigma[i] = 0
        else
            sigma[i] = (f_plus-f_neg)^2/(f_plus+f_neg)
        end
    end

    return sum((sigma.-charge_asym).^2)/length(sigma)
end

### DEPRECATED need later computation of norms
function computeKappa(charges::Vector{R}, charge_asym) where {R<: Real}
    if length(charges)<6
        return -1.0, -1.0
    end 

    value =  zero(R)
    for width in [5,6]
        #min_seq , max_seq = createMinMaxSequences(charges; InsertEmpty=false)
        #if length(min_seq)<6
        #    append!(min_seq, [0,0,0,0,0,0])
        #end

        
        #kappa_min = computeKappa_Helper(min_seq, charge_asym;width=width)
        #kappa_max = computeKappa_Helper(max_seq, charge_asym;width=width)
        kappa = computeKappa_Helper(charges, charge_asym;width=width)

        kappa_norm = kappa#/kappa_max#(kappa .-kappa_min)./(kappa_max-kappa_min) 

        value += kappa_norm
    end
    return value/2.0
end

function computeChargeStuff(charges::Vector{R}) where {R<: Real}
    N = length(charges)
    q_total = sum(charges)

    q_pos = zero(R)
    N_pos = zero(Int32)
    q_neg = zero(R)
    N_neg = zero(Int32)

    for q in charges
        if q >0 
            q_pos += q 
            N_pos += 1
        elseif q<0 
            q_neg += q 
            N_neg += 1
        end
    end

    f_plus =  N_pos/N
    f_neg  =  N_neg/N
    fcr = f_plus+f_neg
    ncpr = abs(f_plus - f_neg)

    charge_asym = (f_plus-f_neg)^2/(f_plus+f_neg)

    return N, q_total, q_pos, q_neg, f_plus, f_neg, fcr, ncpr, charge_asym
end

function AnalyseChargeSequence(charges::Vector{R}) where {R<: Real}


    scd = SequenceChargeDecorator(charges)
    nscd = 0.0 #normSCDValues([scd], charges)[1]

    scd_V1 = SequenceChargeDecorator_Var1(charges)
    #nscd_V1 = normSCDValues([scd_V1], charges; func=SequenceChargeDecorator_Var1)[1] ### not properly normalised

    scd_V2 = SequenceChargeDecorator_Var2(charges)
    #nscd_V2 = normSCDValues([scd_V2], charges; func=SequenceChargeDecorator_Var2)[1] ### not properly normalised

    ### taken from "Conformations of intrinsically disordered proteins are influenced by linear sequence distributions of  oppositely charged residues
    N, q_total, q_pos, q_neg, f_plus, f_neg, fcr, ncpr, charge_asym = computeChargeStuff(charges)
    
    kappa = computeKappa(charges, charge_asym)

    return N , q_total, q_pos, q_neg, f_plus, f_neg, fcr, ncpr, scd, nscd, scd_V1, scd_V2, charge_asym, kappa
end  