# This code was created in collaboration with my colleague Miguel La Iglesia Mirones
# Lucía Muñoz Gil, MSc in Computational Biology, Bioinformatics Programming Challenges course, 2023
# Code for Assigment #1

# Code for defining the "SeedStock" Class Object

class SeedStock

    # Accesor to set/get the values of the properties from outside the Object
    attr_accessor :seed_stock_ID, :mutant_gene_ID, :last_planted, :grams_remaining, :storage
 
    # Initialize method defining the main properties and setting default values when an instance is created
    def initialize(
        seed_stock_ID: "XXXX",
        mutant_gene_ID: "ATXXXXX",
        last_planted: "dd/mm/yyyy",
        storage:"camaX",
        grams_remaining: 0)

        @seed_stock_ID = seed_stock_ID  # Stock ID
        @mutant_gene_ID = mutant_gene_ID    # Gene ID
        @last_planted = last_planted    # Date of last planted
        @storage = storage  # Storage place
        @grams_remaining  = grams_remaining #  Grams of seed remaining after the last planting
    end

    # Method for simulating the planting of the seeds. Recieves as input the quantinty of grams to be planted
    def plant_seeds(gr_planted)
        @last_planted = Date.today.strftime("%e/%m/%Y") # Update the last_planted property with the date of today (in the same format)
        if @grams_remaining > gr_planted # If the quantity remaining is greater than the one to be planted then
            @grams_remaining -= gr_planted  # reduce by -gr_planted the grams_remaining property
        else    # If there are exactly "gr_planted" gr or less then the remaining stock is 0 gr (so only the remaining stock could be planted)
            warn "WARNING: we have run out of Seed Stock #{@seed_stock_ID}. #{@grams_remaining} grams were planted" 
            @grams_remaining = 0 
        end
    end
end