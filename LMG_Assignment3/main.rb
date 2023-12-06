#----------------- LOAD EXTERNAL LIBRARIES AND MODULES-------------------------------------------------------------------
require 'bio'
require 'net/http'
require  './EmblEntry.rb'
require './EmblProcessor.rb'

#------------------------------ INPUT FILES NAMES AS ARGUMENTS FROM  COMMAND LINE ---------------------------------------

# Check for one common error: not specifying the correct number of files needed for the program to run
if ARGV.length != 1
    abort "Incorrect number of files passed. Two files names must be specified: input list of genes and name for final_report"
end
  
# Check for second common error: incorrect usage, files in incorrect order or wrong name passed.
if ARGV[0] == "ArabidopsisSubNetwork_GeneList.txt"
#if ARGV[0] == 'soloungen.txt'
        gene_file = ARGV[0]
else 
    abort "Incorrect order of files passed. Usage: ruby main.rb ArabidopsisSubNetwork_GeneList.txt"
end


# ------------ PUBLIC INSTANCE METHODS ---------------------------------------

# Get response from URL
def fetch(uri_str) 
    address = URI(uri_str)  
    response = Net::HTTP.get_response(address)
    case response
      when Net::HTTPSuccess then
        return response.body
      else
        raise Exception, "Something went wrong... the call to #{uri_str} failed; type #{response.class}"
        response = false
        return response
    end
end

# Get response with NCBI E-Utils, given a database, file format and gene locus
def ncbi_fetch(database, file_format, id)
    url = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=#{database}&format=#{file_format}&id=#{id}"
    #puts url
    res = fetch(url)
    return res
end


# Transform locus names file into string of coma separated locus names, while checking A. thaliana format
def read_from_file(filename)
    gene_list = []
    gene_file = File.open(filename, 'r')
    gene_file.readlines.each do |line|  
        locus_name=line.chomp
        if locus_name !~ /A[Tt]\d[Gg]\d\d\d\d\d/ # Check for correct format of gene locus name
            abort "Locus name #{locus_name} does not meet the correct format. Please define locus names as AT0g00000"
        end
        gene_list << locus_name 
    end
    gene_file.close
    gene_ids = gene_list.join(",")
    return gene_ids
end

# Gets coordinates from string
def get_coordinates(coordinates_string)
    position = coordinates_string.match(/\d+\.\.\d+/)[0]
    start_coordinate = position.split('..')[0].to_i
    end_coordinate = position.split('..')[1].to_i
    return[start_coordinate, end_coordinate]
end

# Given an specific motif (sequence pattern) get start and end coordinates of any match in a given sequence
def find_motifs (sequence, motif)
    matches = []
    position = 0
    while (match = sequence.match(motif, position))
        match_sequence = match[0]
        start_motif = match.begin(0)
        end_motif = match.end(0)
        matches << {
          sequence: match_sequence,
          start: start_motif.to_i,
          end: end_motif.to_i
        }
        # Move the position forward, "- 5" allows to detect a domain "CTTCTTCTT", and that will be detected as 2
        position = end_motif - 5
    end
    return matches
end

# Write a GFF3 file
def write_gff(gff, coordinates_to_use)
    file = coordinates_to_use == 'chromosomal coordinates' ? 'AT_repeats_chromosomal.gff' : 'AT_repeats_sequence.gff'
    File.open(file, 'w') do |file|
      file.puts gff.to_s
    end
end
  
#-------------------------------------------MAIN CODE ------------------------------------------------------------------
# BEFORE EVERYTHING, SEARCHING EMBL FILE (Task 1)

# 1:  Using BioRuby, examine the sequences of the ~167 Arabidopsis genes from the last assignment by retrieving them from whatever database you wish #
# ArabidopsisSubNetwork file -> each gene locus (lines) -> search for embl file -> add to 'AT_sequences.embl'

puts "Processing #{gene_file} file, this might take a while..."

# ACTIVAR AL FINAL
# Create string with gene IDs from file
#gene_ids = read_from_file(gene_file)

# Only one search with list of IDs
#response_body = ncbi_fetch(database = 'ensemblgenomesgene',file_format = 'embl', id = gene_ids)

#output_file = File.open('AT_sequences.embl', 'w')
#output_file.write(response_body)
#output_file.close


search_positive = Bio::Sequence::NA.new("cttctt")
search_complementary = Bio::Sequence::NA.new("aagaag")
REGEX_POSITIVE = Regexp.new(search_positive.to_re)
REGEX_COMPLEMENT = Regexp.new(search_complementary.to_re)


# 2, 3, 4a and 4b

eobject_processor = EmblProcessor.new('AT_sequences.embl')
gff_seq = eobject_processor.load_to_gff('sequences coordinates')
write_gff(gff_seq, 'sequences coordinates')

eobject_processor_chr = EmblProcessor.new('AT_sequences.embl')
gff_chr = eobject_processor_chr.load_to_gff('chromosomal coordinates')
puts gff_chr
write_gff(gff_chr, 'chromosomal coordinates')


