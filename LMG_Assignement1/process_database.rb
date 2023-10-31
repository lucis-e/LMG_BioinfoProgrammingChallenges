# This code was created in collaboration with my colleague Miguel La Iglesia Mirones
# Lucía Muñoz Gil, MSc in Computational Biology, Bioinformatics Programming Challenges course, 2023
# Code for Assigment #1

#--------LOAD EXTERNAL LIBRARIES AND MODULES--------------------------------------------------

require 'date' # Date library to update the date of the last seed planting
require './SeedStock' # Load all 4 classes created for this assignement
require './HybridCross'
require './Gene'
require './SeedStockDatabase'


#----INPUT FILE NAMES AS ARGUMENTS FROM COMMAND LINE-------------------------------------------

if ARGV.length != 4 # Check for one common human error: not specifying all 4 file names needed for the assignement
    puts "Incorrect number of files passed. There must be 4 files"
    exit    # if the number of files is not 4, then the program stops running
end

# Check for another common human error: not specifying the files in the correct order. If all files are in the correct
# position then declare specific variables and the program keeps runing. If the order is wrong, the script stops.
if ARGV[0] == "gene_information.tsv" && ARGV[1] == "seed_stock_data.tsv" && ARGV[2] == "cross_data.tsv" && ARGV[3] == "new_stock_file.tsv"
    gene_info = ARGV[0]
    seed_stock_file = ARGV[1]
    cross_data = ARGV[2]
    new_file_path = ARGV[3]
else 
    puts "Incorrect order of files passed. Usage: ruby process_database.rb gene_information.tsv seed_stock_data.tsv cross_data.tsv new_stock_file.tsv"
    exit
end


#---------------- MAIN CODE--------------------------------------------------------------------

# TASK 1: Seed planting simulation and new genebank state printed to new file

sesion = SeedStockDatabase.new  # Create an instance of the class "SeedStockDatabase". This instance would allow reading
                                # and writing files, as well as keeping record of the "SeedStock" and "Gene" class 
                                # instances created.

sesion.load_from_file(seed_stock_file)  # Load info from seed_stock_data.tsv file and create new SeedStock objects, which
                                        # are saved into hash (and would be accesible by SeedStock ID). Also header
                                        # of file is saved (for using it when writing new file).

sesion.stock_instances.each do |_key, gene_object|  # For each SeedStock object in previously created hash, the
    gene_object.plant_seeds(7)                      # planting of 7 gr of seeds is simulated and genebank info is
end                                                 # updated. Hash has the updated info as it storages class objects

sesion.write_database(new_file_path)  # New genbank state (represented in the SeedStock objects) is printed into new file

# TASK 2: Chi-squared test for searching for linked genes and adding it as new gene property

sesion.load_from_file(gene_info)    # Load info from gene_information.tsv file and create new Gene objects, which 
                                    # are saved into hash (and are accesible by Gene ID).

File.open(cross_data, "r").each.with_index do |line, line_num|  # read the information in cross_data.tsv file line by line
    next if line_num==0 # skip header

    # Create new HybridCross instance
    p1, p2, f2_w, f2p1, f2p2, f2p1p2 = line.split() # each line is splited into array ( by "\t" separator) and each element
                                                    # of the array is assigned to a variable

    crossing = HybridCross.new(parent1_ID:p1,   # Create new HybridCross instance (one per line) and define properties 
        parent2_ID:p2,                          # based on the values retrieved from reading the specific line
        f2_wildtype:f2_w.to_f,
        f2_p1:f2p1.to_f,
        f2_p2:f2p2.to_f,
        f2_p1p2:f2p1p2.to_f)

    x_squared = crossing.chi_squared()  # Perform the chi_squared test with Chi_squared method (HybridCross method).
                                        # This method would calculate and return the chi-square value based on the F2 data

    if x_squared > 7.82 # Critical chi-square value for a test of significance at alpha = 0.05 and 3 degrees of freedom
                        # 3 d.o.f. since there are 4 grous (wildtipe, mutated gen 1, mutated gen 2, both mutated)
                        # Null hypothesis is genes not being linked (the F2 have equal probability of inheriting all 
                        # possible genotypic combinations). H0 can be rejected if chi-square value is greater 7.82.

        # If H0 is rejected (x_squared > 7.82), then genes are linked. To get the name of the genes that are linked:
        # 1. Get the SeedStock objects with stock ID that matches each of the parents IDs and extract the gene IDs
        linked_gene1 = sesion.get_seed_stock(crossing.parent1_ID).mutant_gene_ID
        linked_gene2 = sesion.get_seed_stock(crossing.parent2_ID).mutant_gene_ID

        # 2, Get the Gene objects that match the gene IDs of the parents so we have all the information of the gene
        linked_genename1 = sesion.get_gene_info(linked_gene1)
        linked_genename2 = sesion.get_gene_info(linked_gene2)

        # 3. Add as property of the Gene object (values are the names of genes that each gene is linked to)
        linked_genename1.linked_to = linked_genename2.name # add new value (linked gene name) to the property
        linked_genename2.linked_to = linked_genename1.name
        
        puts "Recording: #{linked_genename1.name} is genetically linked to #{linked_genename2.name} with chisquare score #{x_squared}"
    end
end

# Finally, a final report of linked genes discovery is displayed:
puts
puts
puts "Final report: \n"
puts

sesion.gene_instances.each do |_key, gene_object|   # for each Gene object storaged in the hash:
    if gene_object.respond_to?(:linked_to) && gene_object.linked_to # check if the object has the property linked_to
                                                                    # (if the property has values asigned then the gene
                                                                    # is linked to some other gene(s))
        puts "#{gene_object.name} is linked to #{gene_object.linked_to}"    # Display only linkage info of those genes
                                                                            # that we discovered that are actually linked
                                                                            # to some other gene(s)
    end
end
