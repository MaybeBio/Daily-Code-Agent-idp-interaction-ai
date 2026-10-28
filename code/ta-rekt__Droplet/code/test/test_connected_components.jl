# MARK: -CC multi
    @testset "Connected components (multiple)" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = BulkMonomer(;
            length = 2,
            initShape = "box",
            initParams = 2
        )

        fb = FastBulkMonomer(;
            length = 2,
            initShape = "box",
            initParams = 2
        )

        # Set interaction potentials after init
        intPot!(b, -1.0)
        intPot!(fb, -1.0)

        init!(l, [b, fb])

        @test cc(b) == Vector([Set([4,3]), Set([2,1])])
        @test cc(fb) == Vector([Set([4,3]), Set([2,1])])
    end

# MARK: -CC single
    @testset "Connected components (single)" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = BulkMonomer(;
            length = 20,
            initShape = "box",
            initParams = 1
        )

        fb = FastBulkMonomer(;
            length = 20,
            initShape = "box",
            initParams = 1
        )

        # Set interaction potentials after init
        intPot!(b, -1.0)
        intPot!(fb, -1.0)

        init!(l, [b, fb])

        fx = cc(fb)
        fy = Vector{Set{Int32}}()
        fy0 = Set{Int32}(range(1,20))
        push!(fy, fy0)
    
        @test fx == fy

        x = cc(b)
        y = Vector{Set{Int32}}()
        y0 = Set{Int32}(range(1,20))
        push!(y, y0)
    
        @test x == y
    end

# MARK: -CC line
    @testset "Connected components (line)" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 2
        )

        b = BulkMonomer(;
            length = 8,
            initShape = "line"
        )

        fb = FastBulkMonomer(;
            length = 8,
            initShape = "line"
        )

        # Set interaction potentials 
        intPot!(b, -1.0)
        intPot!(fb, -1.0)

        init!(l, [b, fb])

        x = cc(b)
        y = Vector{Set{Int32}}()
        y0 = Set{Int32}(range(1,8))
        push!(y, y0)

        @test x == y

        fx = cc(fb)
        fy = Vector{Set{Int32}}()
        fy0 = Set{Int32}(range(1,8))
        push!(fy, fy0)

        @test fx == fy
    end

