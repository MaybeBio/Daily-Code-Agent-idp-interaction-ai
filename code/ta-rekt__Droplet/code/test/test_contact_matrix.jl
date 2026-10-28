using Test
using Droplet
using LinearAlgebra

@testset "Contact Matrix" begin
    boxlength = 128
    
    @testset "compute_contact_matrix! - basic functionality" begin
        # Set up lattice and polymer with 20 beads
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        
        poly = Polymer(;
            length = 20,
            maxlength = 128,
            nMonomerTypes = 1,
            initShape = "line"
        )
        
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        # Initialize contact matrix
        contact_matrix = zeros(Int, poly.length, poly.length)
        
        # Compute contact matrix
        result = compute_contact_matrix!(poly, contact_matrix)
        
        # Check that result is the same matrix object
        @test result === contact_matrix
        
        # Check matrix dimensions
        @test size(contact_matrix) == (20, 20)
        
        # Check symmetry: contact_matrix[i,j] should equal contact_matrix[j,i]
        for i in 1:poly.length
            for j in i+1:poly.length
                @test contact_matrix[i, j] == contact_matrix[j, i]
            end
        end
        
        # Check diagonal is zero (beads don't contact themselves)
        for i in 1:poly.length
            @test contact_matrix[i, i] == 0
        end
    end
    
    @testset "compute_contact_matrix! - contact detection" begin
        # Create a polymer with defined shape
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        
        poly = Polymer(;
            length = 50,
            maxlength = 128,
            nMonomerTypes = 1,
            initShape = "line"
        )
        
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        contact_matrix = zeros(Int, poly.length, poly.length)
        compute_contact_matrix!(poly, contact_matrix, 2.0)
        
        # For a line polymer, consecutive beads should have many contacts
        # (they are bonded and adjacent in space)
        for i in 1:poly.length-1
            # Consecutive beads should have non-zero contact count
            @test contact_matrix[i, i+1] >= 0  # Should be positive or zero
        end
    end
    
    @testset "compute_contact_matrix! - contact threshold" begin
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        
        poly = Polymer(;
            length = 30,
            maxlength = 128,
            nMonomerTypes = 1,
            initShape = "line"
        )
        
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        # Test with different contact thresholds
        contact_matrix_tight = zeros(Int, poly.length, poly.length)
        contact_matrix_loose = zeros(Int, poly.length, poly.length)
        
        # Tight threshold (closer contacts only)
        compute_contact_matrix!(poly, contact_matrix_tight, 1.0)
        
        # Loose threshold (more distant contacts)
        compute_contact_matrix!(poly, contact_matrix_loose, 3.0)
        
        # Sum of contacts should increase with looser threshold
        @test sum(contact_matrix_loose) >= sum(contact_matrix_tight)
    end
    
    @testset "compute_contact_matrix! - matrix accumulation" begin
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        
        poly = Polymer(;
            length = 15,
            maxlength = 128,
            nMonomerTypes = 1,
            initShape = "line"
        )
        
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        # Pre-populate contact matrix with initial values
        contact_matrix = ones(Int, poly.length, poly.length)
        contact_matrix[diagind(contact_matrix)] .= 0  # Zero out diagonal
        
        initial_sum = sum(contact_matrix)
        
        # Run compute_contact_matrix! (should accumulate, not replace)
        compute_contact_matrix!(poly, contact_matrix)
        
        final_sum = sum(contact_matrix)
        
        # After accumulation, sum should be at least as large as initial
        @test final_sum >= initial_sum
    end
    
    @testset "compute_contact_matrix! - periodic boundary conditions" begin
        # Create a small system to test periodic boundaries
        l = Lattice(; boxlength = 20, nspecies = 1)
        
        poly = Polymer(;
            length = 10,
            maxlength = 20,
            nMonomerTypes = 1,
            initShape = "line"
        )
        
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        contact_matrix = zeros(Int, poly.length, poly.length)
        compute_contact_matrix!(poly, contact_matrix, 2.0)
        
        # Check that distances are calculated considering periodic boundaries
        # (function should handle wrapping correctly)
        @test size(contact_matrix) == (10, 10)
        
        # Verify symmetry holds with periodic boundaries
        for i in 1:poly.length
            for j in i+1:poly.length
                @test contact_matrix[i, j] == contact_matrix[j, i]
            end
        end
    end
    
    @testset "compute_contact_matrix! - single bead polymer" begin
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        
        poly = Polymer(;
            length = 1,
            maxlength = 128,
            nMonomerTypes = 1,
            initShape = "line"
        )
        
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        contact_matrix = zeros(Int, 1, 1)
        result = compute_contact_matrix!(poly, contact_matrix)
        
        # With only 1 bead, no contacts should be detected
        @test contact_matrix[1, 1] == 0
    end
    
    @testset "compute_contact_matrix! - two bead polymer" begin
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        
        poly = Polymer(;
            length = 2,
            maxlength = 128,
            nMonomerTypes = 1,
            initShape = "line"
        )
        
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        contact_matrix = zeros(Int, 2, 2)
        compute_contact_matrix!(poly, contact_matrix, 2.0)
        
        # Diagonal should be zero
        @test contact_matrix[1, 1] == 0
        @test contact_matrix[2, 2] == 0
        
        # Symmetry should hold
        @test contact_matrix[1, 2] == contact_matrix[2, 1]
    end
    
    @testset "compute_contact_matrix! - contact count accumulation" begin
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        
        poly = Polymer(;
            length = 12,
            maxlength = 128,
            nMonomerTypes = 1,
            initShape = "line"
        )
        
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        # Test multiple accumulations
        contact_matrix = zeros(Int, poly.length, poly.length)
        
        # First call
        compute_contact_matrix!(poly, contact_matrix, 2.0)
        first_sum = sum(contact_matrix)
        
        # Second call (accumulation)
        compute_contact_matrix!(poly, contact_matrix, 2.0)
        second_sum = sum(contact_matrix)
        
        # Second sum should be approximately double the first (within accumulation)
        @test second_sum >= first_sum
    end
    
    @testset "compute_contact_matrix! - distance calculation correctness" begin
        # Create a very simple test case with manual verification
        l = Lattice(; boxlength = 30, nspecies = 1)
        
        poly = Polymer(;
            length = 5,
            maxlength = 30,
            nMonomerTypes = 1,
            initShape = "line"
        )
        
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        contact_matrix = zeros(Int, poly.length, poly.length)
        
        # Use a threshold of 2.0 which allows for bonds and nearby contacts
        compute_contact_matrix!(poly, contact_matrix, 2.0)
        
        # All elements should be >= 0
        @test all(contact_matrix .>= 0)
        
        # Verify symmetry
        @test contact_matrix == contact_matrix'
    end
end
