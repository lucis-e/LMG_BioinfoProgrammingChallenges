#----------------- LOAD EXTERNAL LIBRARIES AND MODULES-------------------------------------------------------------------
require 'bio'
require 'net/http'


#------------------------------ INPUT FILES NAMES AS ARGUMENTS FROM  COMMAND LINE ---------------------------------------

# Check for one common error: not specifying the correct number of files needed for the program to run
if ARGV.length != 1
    abort "Incorrect number of files passed. Two files names must be specified: input list of genes and name for final_report"
end
  
# Check for second common error: incorrect usage, files in incorrect order or wrong name passed.
if ARGV[0] == "ArabidopsisSubNetwork_GeneList.txt"
    input_gene_list = ARGV[0]
else 
    abort "Incorrect order of files passed. Usage: ruby main.rb ArabidopsisSubNetwork_GeneList.txt"
end

#------------------------ PUBLIC INSTANCE METHODS -----------------------------------------------------------------------

# esto lo podemos poner dentro de la de abajo me da igual
def get_sequences_from_ensembl(locus_name)
    address = URI("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{locus_name}")  
    response = Net::HTTP.get_response(address)  # use the Net::HTTP object "get_response" method
    return(response.body)
end

# read the gene locus names from file and retrieve the sequence info from Ensemble, save all records in a local file
def read_from_file(filename)
    gene_file = File.open(filename, 'r')
    gene_file.readlines.each do |line|
        locus_name=line.chomp

        if locus_name !~ /A[Tt]\d[Gg]\d\d\d\d\d/ # Check for correct format of gene locus name
            abort "Locus name #{locus_name} does not meet the correct format. Please define locus names as AT0g00000"
        end

        record = get_sequences_from_ensembl(locus_name)
        File.open('AT_sequences.embl', 'w') do |myfile|  # w makes it writable
            myfile.puts record  # save multiple records in one same file
        end
    end
end

#-------------------------------------------MAIN CODE ------------------------------------------------------------------
puts "Processing #{input_gene_list} file, this might take a while..."

# retrive sequence from Ensembl and save into local file "AT_sequences.embl"

read_from_file(input_gene_list)