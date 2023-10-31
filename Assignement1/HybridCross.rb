# This code was created in collaboration with my colleague Miguel La Iglesia Mirones
# Lucía Muñoz Gil, MSc in Computational Biology, Bioinformatics Programming Challenges course, 2023
# Code for Assigment #1

# Code for defining the "HybridCross" Class Object

class HybridCross

    # Accesor to set/get the values of the properties from outside the Object
    attr_accessor :parent1_ID, :parent2_ID, :f2_wildtype, :f2_p1, :f2_p2, :f2_p1p2

    # Initialize method defining the main properties and setting default values when an instance is created
    def initialize( 
        parent1_ID: "XXXX",
        parent2_ID: "XXXX",
        f2_wildtype: 0,
        f2_p1: 0,
        f2_p2: 0,
        f2_p1p2: 0)

        @parent1_ID = parent1_ID    # Stock ID of the first parent (recessive homocygote mutant for gene specified in SeedStock instance)
        @parent2_ID = parent2_ID    # Stock ID of the second parent (recessive homocygote mutant for gene specified in SeedStock instance)
        @f2_wildtype = f2_wildtype  # Number of F2 individuals displaying wildtype phenotype
        @f2_p1 = f2_p1  # Number of F2 individuals displaying mutant phenotype for same gene as parent 1
        @f2_p2 = f2_p2 # Number of F2 individuals displaying mutant phenotype for same gene as parent 2
        @f2_p1p2  = f2_p1p2 # Number of F2 individuals displaying mutant phenotype for both genes
    end


    # Method for calculating the chi-squared with the F2 data
    def chi_squared()
        total_seeds = @f2_wildtype + @f2_p1 + @f2_p2 + @f2_p1p2 # Total number of F2 individuals 
        
        # Calculate the expected values (multiplying by the expected proportion when 2 genes are not linked)
        exp_w = total_seeds * 9/16 # 9 out of 16 individuals are expected to present wildtype phenotype 
        exp_het_mutant = total_seeds * 3/16 # 3 out of 16 individuals are expected to present mutant phenotype for only one gene and other 3 for the other gene
        exp_hom_mutant = total_seeds * 1/16 # Only 1 individual out of 16 is expected to display double mutant phenotype

        # CalCulate the chi-square
        x_squared = (((@f2_wildtype - exp_w)**2 / exp_w) + ((@f2_p1 - exp_het_mutant)**2 / exp_het_mutant) + ((@f2_p2 - exp_het_mutant)**2 / exp_het_mutant) + ((@f2_p1p2 - exp_hom_mutant)**2 / exp_hom_mutant))
        return x_squared # Return as output the chi-square value of the F2 data
    end
end