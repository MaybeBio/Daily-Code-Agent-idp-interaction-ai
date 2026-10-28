# MARK: -Density uniform
    @testset "Box density (uniform)" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = BulkMonomer(;
            length = 256,
            initShape = "uniform",
            initParams = 2
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        @test boxDensity(b, [[Int32(1), b.boxlength], 
                            [Int32(1), b.boxlength], 
                            [Int32(1), b.boxlength]]) == b.length / b.boxlength ^ 3
    end

# MARK: -Density full
    @testset "Box density (full box)" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = BulkMonomer(;
            length = 27,
            initShape = "box",
            initParams = 1,
            initPos = [2,2,2]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        y = b.length / 3 ^ 3
        x = boxDensity(b, [[Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)]])
        @test x == y

        b = BulkMonomer(;
            length = 26,
            initShape = "box",
            initParams = 1,
            initPos = [2,2,2]
        )

        slist = [b]

        # Set interaction potential
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        y = 26/27
        x = boxDensity(b, [[Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)]])
        @test x == y

        y = 7/27
        x = boxDensity(b, [[Int32(2), Int32(4)], 
                        [Int32(2), Int32(4)], 
                        [Int32(2), Int32(4)]])

        @test x == y

        y = 0.0
        x = boxDensity(b, [[Int32(3), Int32(5)], 
                        [Int32(3), Int32(5)], 
                        [Int32(3), Int32(5)]])

        @test x == y

        y = 0.0
        x = boxDensity(b, [[Int32(4), Int32(6)], 
                        [Int32(4), Int32(6)], 
                        [Int32(4), Int32(6)]])

        @test x == y


        b = BulkMonomer(;
            length = 27,
            initShape = "box",
            initParams = 1,
            initPos = [2,2,1]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        y = 1.0
        x = boxDensity(b, [[Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)], 
                        [Int32(0), Int32(2)]])
        @test x == y

        y = 2/3
        x = boxDensity(b, [[Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)], 
                        [Int32(1), Int32(3)]])
        @test x == y
    end

# MARK: -Split droplet
    @testset "Split droplet" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = BulkMonomer(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,2]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)
    
        ccs = cc(b)
        dlist = splitDroplet!(b, ccs[1])
        @test dlist == [Set([5, 6]), Set([4, 3]), Set([7, 8]), Set([2, 1])]

        b = BulkMonomer(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)
    
        # initialize lattice and species
        init!(l, slist)
    
        ccs = cc(b)
        dlist = splitDroplet!(b, ccs[1])
        @test dlist == [Set([5]), Set([4]), Set([6]), Set([7]), 
                      Set([2]), Set([8]), Set([3]), Set([1])]

        b = BulkMonomer(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [1,2,2]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)
    
        ccs = cc(b)
        dlist = splitDroplet!(b, ccs[1])
        
        @test dlist == [Set([5, 6, 7, 8]), Set([4, 2, 3, 1])]

        b = BulkMonomer(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [2,2,2]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)
    
        ccs = cc(b)
        dlist = splitDroplet!(b, ccs[1])

        @test dlist == [Set([5, 4, 6, 7, 2, 8, 3, 1])]
    end

# MARK: -Droplet bbox
    @testset "Droplet bbox" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = BulkMonomer(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [8,8,8]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)
        bbox = dropletBbox(b, ccs[1])

        @test bbox == [[6, 9], [6, 9], [6, 9]]

        b = BulkMonomer(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [2,2,2]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)
        bbox = dropletBbox(b, ccs[1])

        @test bbox == [[1,3],[1,3],[1,3]]

        b = BulkMonomer(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [3,2,2]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)
        bbox = dropletBbox(b, ccs[1])

        @test bbox == [[1,4],[1,3],[1,3]]
    end

# MARK: -Fill/clear
    @testset "Fill/clear bbox" begin

        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        for i in range(1, 100)

            b = BulkMonomer(;
                length = i,
                initShape = "box",
                initParams = 1,
                initPos = [3,3,3]
            )

            slist = [b]

            # Set interaction potential 
            intPot!(b, -1.0)

            # initialize lattice and species
            init!(l, slist)

            ccs = cc(b)

            bbox = dropletBbox(b, ccs[1])

            fillBbox!(b, ccs[1], bbox)

            @test length(collect(keys(b.ghost.pos))) == vol(bbox) - length(ccs[1])
            @test sum(map(x -> x > 0, b.ghost.latS)) == vol(bbox) - length(ccs[1])

            clearBbox!(b, bbox)

            @test iszero(b.ghost.latS)
            @test isempty(b.ghost.pos)
        end
    end
# MARK: -Droplet volume
    @testset "Droplet volume" begin
        l = Lattice(;
            boxlength = 16,
            nspecies = 1
        )

        b = BulkMonomer(;
            length = 64,
            initShape = "box",
            initParams = 1,
            initPos = [3,3,3]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        keys_ = [[2,2,2],[2,2,3],[2,3,2],[2,3,3],
             [3,2,2],[3,2,3],[3,3,2],[3,3,3],[1,1,1],[2,1,2]]  # poke a hole in it
        for key in keys_
            x, y, z = key
            key = b.latS[x,y,z]
            b.latS[x,y,z] = 0
            pop!(b.pos, key)
        end

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == 54

        b = BulkMonomer(;
            length = 64,
            initShape = "box",
            initParams = 1,
            initPos = [3,3,3]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        keys_ = [[2,2,2],[2,2,3],[2,3,2],[2,3,3],
             [3,2,2],[3,2,3],[3,3,2],[3,3,3],[1,1,1]]  # no hole: outer volume preserved
        for key in keys_
            x, y, z = key
            key = b.latS[x,y,z]
            b.latS[x,y,z] = 0
            pop!(b.pos, key)
        end

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == 63


        b = BulkMonomer(;
            length = 64,
            initShape = "box",
            initParams = 1,
            initPos = [3,3,3]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        keys_ = [[2,2,2],[2,2,3],[2,3,2],[2,3,3],
                 [3,2,2],[3,2,3],[3,3,2],[3,3,3]]
        for key in keys_
            x, y, z = key
            key = b.latS[x,y,z]
            b.latS[x,y,z] = 0
            pop!(b.pos, key)
        end

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == 64

        b = BulkMonomer(;
            length = 24,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == length(collect(keys(b.pos)))


        b = BulkMonomer(;
            length = 29,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == length(collect(keys(b.pos)))


        b = BulkMonomer(;
            length = 2,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == length(collect(keys(b.pos)))

        b = BulkMonomer(;
            length = 8,
            initShape = "box",
            initParams = 1,
            initPos = [1,1,1]
        )

        slist = [b]

        # Set interaction potential 
        intPot!(b, -1.0)

        # initialize lattice and species
        init!(l, slist)

        ccs = cc(b)

        @test dropletVolume(b, ccs[1]) == length(collect(keys(b.pos)))
    end
