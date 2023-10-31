# This code was created in collaboration with my colleague Miguel La Iglesia Mirones
# Lucía Muñoz Gil, MSc in Computational Biology, Bioinformatics Programming Challenges course, 2023
# Code for Assigment #1

# Code for defining the "SeedStockDatabase" Class Object

class SeedStockDatabase

    # Accesor to set/get the values of the properties from outside the Object
    attr_accessor :stock_instances, :gene_instances

    # Initilize method defining the main properties
    def initialize
        @stock_instances = {}   # HASH that allows us to store and access (by Stock ID) all SeedStock instances created
        @gene_instances = {}    # Same as before but to store and access (by Gene ID) all Gene instances crated
        @header = "NA" # Will store the header of the file readed by load_from_file method and used in write_database method
    end

    # Method for indexing SeedStock instances into the stock_instances hash
    def add_seed_stock(seed)
        @stock_instances[seed.seed_stock_ID] = seed # Stock ID is the key and SeedStock object the value
    end

    # Method for retrieving SeedStock objects by Stock ID
    def get_seed_stock(seed_stock_ID)
        return @stock_instances[seed_stock_ID]  # Return the object (value) for key "seed_stock_ID"
    end

    # Method for indexing Gene instances into the gene_instances has
    def add_gene(gene)
        @gene_instances[gene.gene_ID] = gene    # Gene ID is the key and Gene object the value
    end

    # Method for retrieving Gene objects by Gene ID
    def get_gene_info(gene_ID)
        return @gene_instances[gene_ID] # Return the object (value) for key "gene_ID"
    end

    # Method for loading data from files
    def load_from_file(file)

        # This code would act different depending on which file it is passed as argument (want to be loaded)
        File.open(file, "r").each.with_index do |line, line_num|    # with_index will alow us to filter the first row (header) by row index (index = 0)
                                                                    # .each allows us to read file content line by line
        if file == "seed_stock_data.tsv"    # If the file to be loaded is "seed_stock_data.tsv" then:
            if line_num==0 # save header in "header" variable and skip to next file line
                @header = line
                next
            end

            seed_ID, gene_id, last_plant_date, storage_place, grams = line.split()  # split file line (by "\t" separator) into array 
                                                                                    # and each element is assigned to a variable

            seed_record = SeedStock.new(    # Create new SeedStock instance (one per line) and set values for its properties based on
                seed_stock_ID:seed_ID,      # the values retrieved from reading that specific line
                mutant_gene_ID:gene_id,
                last_planted:last_plant_date,
                storage:storage_place,
                grams_remaining:grams.to_i)
            
            add_seed_stock(seed_record)     # Index newly created SeedStock instance in hash (to search by Stock ID)

        elsif file == "gene_information.tsv" # On the other hand, if the files to be oaded is "gene_information,tsv" then:

            next if line_num==0 # skip header (no need to be storaged since we only want to print the state of the database)
            
            geneid, gene_name, phenotype = line.split() # split file line (by "\t" separator) into array and assing each element to a variable

            gene_record = Gene.new( # Create new Gene instance (one per line) and set values for "name" and "phenotype" properties
                name:gene_name,     # based on the info retrieved from reading that line
                mutant_phenotype:phenotype)
            
            gene_record.gene_ID=(geneid)    # Value for gene_ID property should be setted using the "set" gene function so our program test Gene ID format
            add_gene(gene_record)   #  Index newly created Gene instance in hash (to search by Gene ID)
        end
    end
    end

    # Method for writting the new Database state into new file
    def write_database(newfile) # As input recieves the name of the new file to be created
        File.open(newfile, "w") do |file|   # "w" would create if the file does not exist or re-write if it does
            file.puts @header   # Header retrieved from loading the seed_stock_data.tsv is written as first line
            @stock_instances.each_value do |seed|   # For each value of the hash storing the SeedStock instances, print each of its properties in the same line
                file.puts "#{seed.seed_stock_ID}\t#{seed.mutant_gene_ID}\t#{seed.last_planted}\t#{seed.storage}\t#{seed.grams_remaining}"
            end
        end
    end
end