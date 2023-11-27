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
        File.open('AT_sequences.embl', 'a') do |myfile|  # a open the file in append mode, create if not exist, append if exists
            myfile.puts record  # save multiple records in one same file
        end
    end
end

#-------------------------------------------MAIN CODE ------------------------------------------------------------------
puts "Processing #{input_gene_list} file, this might take a while..."

# 1. retrive sequence from Ensembl and save into local file "AT_sequences.embl"

#read_from_file(input_gene_list) esto solo lo activamos la ult vez que ya tenemos el archivo creado

# 2. Loop over very exon feature and scan it for the CTTCTT sequence
datafile = Bio::FlatFile.auto('AT_sequences.embl')
puts datafile.class  # Bio::FlatFile

#datafile.each_entry do |entry| # the FILE is not the same as the RECORD - multiple records can exist in a file
#    next unless entry.accession     # scape empty and nil values
#  
#end

entry=datafile.next_entry

puts entry.class
puts "# #{entry.accession} - #{entry.species}"
entry.features.each do |feature|
    position = feature.position
    qual = feature.assoc            # feature.assoc gives you a hash of Bio::Feature::Qualifier objects 
                                    # i.e. qualifier['key'] = value  for example qualifier['gene'] = "CYP450")

    next unless qual['exon']    # this is an indication that the feature is a transcript

    # collects gene name and exon positions and joins it into a string
    gene_info = [qual['gene'], qual['exon']].compact.join(', ')
    puts gene_info
end

