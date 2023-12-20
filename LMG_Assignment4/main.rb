require 'bio'

# -----------------------CONSTANTS ----------------------- #

# Start codon of a CDS sequence
START_CODON = 'ATG'
# Stop codon of a CDS sequence
STOP_CODON = ['TAA', 'TAG', 'TGA']
# Evalue threshold for blast report hits
EVALUE_THRESHOLD = 1e-6  

# ----------------------- MAIN FUNCTIONS ----------------------- #

# Check content of the fasta files of the proteomes, the sequences are genomic sequences or protein sequences
# @param fasta [String] fasta file
# @return [String] 'nucl' for nucleotide and 'prot' for protein

def fasta_sequence(fasta)
    first_seq = fasta.first.seq
    biosequence = Bio::Sequence.new(first_seq).auto   # guess the sequence type with .auto
    if biosequence.instance_of? Bio::Sequence::NA
        return 'nucl'
    elsif biosequence.instance_of? Bio::Sequence::AA
        return 'prot'
    end
end 

# Filter out those sequences that do not start or end with the start and one of the stop codons in eukaryotes
# @param nucleotide_fasta_file [String]  fasta file
# @param start_codon [String] consensus start codon
# @param stop_codons [Array] consensus stop codons
# @return [Array] filtered out entries that meet the requirements

def only_cds?(nucleotide_fasta_file, start_codon, stop_codons)
    filtered_proteome = []
    nucleotide_fasta_file.each do |entry|
        next unless entry.seq.start_with?(start_codon)
        next unless stop_codons.any? {|stop_codon| entry.seq.end_with?(stop_codon)}
        filtered_proteome << entry
    end
    return filtered_proteome
end

# Translate CDS to protein: to perform a reciprocal blastp (protein query and protein db)
# @param cds_filtered [Array] filtered array with all entries that start and end with the consensus codons
# @param translated_filename [String] file name of the new translated .fa file
# @return [Bio::FlatFile] fasta file with the translated proteome 

def translate_cds_to_protein(cds_filtered, translated_filename)
    proteome_file = File.open(translated_filename, 'w')
    cds_filtered.each do |entry|
        cds_seq = Bio::Sequence.auto(entry.seq)
        protein_seq = Bio::Sequence.auto(cds_seq.translate.chomp("*"))
        proteome_file.puts protein_seq.output_fasta(entry.definition) #entry.entry_id cuando busque
    end

    proteome_file = Bio::FlatFile.auto(translated_filename)
    return proteome_file
end

# Creates BioRuby objects from fasta files of the proteomes. Filters cds proteome keeping those sequences that start and end with the consensus codons and translates proteome into protein sequences
# @param proteome_file [String] fasta filename to manage
# @return [proteome_file] new fasta file name (for S.pombe is the same, for A.thaliana is the translated file)
# @return [proteome_fasta] BioRuby object containing all sequences

def create_and_check_fasta(proteome_file)
    proteome_fasta = Bio::FlatFile.open(Bio::FastaFormat, proteome_file) 
    sequence_type = fasta_sequence(proteome_fasta)
    case sequence_type
    when 'nucl'
        proteome_fasta = only_cds?(proteome_fasta, START_CODON, STOP_CODON)
        puts "Translating #{proteome_file}. This migth take a while..."
        proteome_file = "TAIR10_translated.fa"
        proteome_fasta = translate_cds_to_protein(proteome_fasta, proteome_file)
    end
    return [proteome_file, proteome_fasta]
end


# Make database from proteome file depending of the type of sequences contained in the fasta file 
# @param fasta_file [Bio::FastaFormat] BioRuby object containing all sequences
# @param filename [String] file name of the proteome to create the database from
# @param dbname [String] DataBase name 
# @return [void]

def make_db(fasta_file, filename, dbname)
    sequence_type = fasta_sequence(fasta_file)    # get the sequence type to define -dbtype
    result = system("makeblastdb -in #{filename} -dbtype '#{sequence_type}' -out ./Databases/#{dbname}")
    if $?.success?
        puts "BLAST database #{dbname} created successfully.\n"
    else
        puts "Error creating BLAST database. Output:\n#{result}\n"
    end
end

# Performs Blast and retrieves the best hit based on e-value
# @param entry [Bio::FastaFormat] the query sequence
# @param factory [Bio::Blast] the database to perform blast against to 
# @return [Bio::Blast] best hit based on e-value

def blast_and_best_hit(entry, factory)
    query = ">#{entry.entry_id}\n#{entry.seq}"
    report = factory.query(query)
  
    filtered_hits = report.hits.select{|hit| !hit.evalue.nil? && hit.evalue <= EVALUE_THRESHOLD} # Blast parameters: maximum e-value threshold 
    sorted_hits = filtered_hits.sort_by {|hit| hit.evalue} # retrieve only the hit with the lowest value

    best_hit = sorted_hits.first
    return best_hit
end

# Search for the original sequence in the proteome to perform second blast (and not use the result of the alignment) given an identifier
# @param identifier [String] protein ID
# @param fasta_file [String] proteome to search the protein
# @return [Bio::FastaForma] complete entry of the protein with the given ID.

def search_sequence(identifier,fasta_file)
    fasta = Bio::FlatFile.open(Bio::FastaFormat, fasta_file)

    found_entry = fasta.find { |entry| entry.entry_id.match?(Regexp.escape(identifier))}    # search for the entry given an ID
    if found_entry
        return found_entry
    else
        puts "COULD NOT FIND YOUR SEQUENCE"
        return
    end
end

# Creates a report of the found putative orthologue candidates
# @param cadidate_hash [Hash] pair of putative orthologue candidates founf
# @param output_report_file [String] Report's file name
# @return [void]

def write_candiates_report(candidate_hash, output_report_file)

    File.open(output_report_file, 'w') do |file|

    file.puts "This code was created by Miguel La Iglesia Mirones and Lucía Muñoz Gil"
    file.puts "Bioinformatics Programming Challenges Course at MSc Computational Biology (UPM)"
    file.puts "December 2023"
    file.puts "-----------------------------------------------------------------------------------------------."
    file.puts
    file.puts "PUTATIVE ORTHOLOGUE CANDIDATES AMONG ARABIDOPSIS AND S.POMBE BY PERFORMING RECIPROCAL BEST BLAST" 
    file.puts
    file.puts "------------------------------------------------------------------------------------------------"
    file.puts
    file.puts "Total number of Orthologue candidate pairs: #{candidate_hash.size}"
    file.puts 
    file.puts "Total number of Orthologue candidates: #{candidate_hash.size * 2}"
    file.puts
    candidate_hash.each do |key, value|
        file.puts "Arabidopsis protein #{key} and S.pombe protein #{value} are Orthologue candidates"
    end
    file.puts
    file.puts "------------------------------------------------------------------------------------------------"
    file.puts
    end
end 

# ------------------------ MAIN CODE --------------------------------------------

# Check input arguments: introducing the file names of the proteomes of the species to discover putative orthologues 

if ARGV.length != 2
    abort "Incorrect number of files passed. Proteome files of the two species should be specified"
end

# Check order of arguments: to make the databases, we should know which file correspond to which species

if ARGV[0] == "TAIR10_cds.fa" && ARGV[1] == "proteome_pombe.fa"    
    arabidopsis_proteome = ARGV[0]
    pombe_proteome = ARGV[1]
elsif ARGV[1] == "TAIR10_cds.fa" && ARGV[0] == "proteome_pombe.fa"
    pombe_proteome = ARGV[0]
    arabidopsis_proteome = ARGV[1]
else
    abort "Incorrect files passed. Usage: ruby main.rb TAIR10_cds.fa proteome_pombe.fa"
end


# 1st: Create BioRuby objects from proteome files and check the ones containing CDS nucleotide sequences 
# Translate nucleotide sequences proteomes into protein sequences proteomes if any

arabidopsis_proteome, arabidopsis_fasta = create_and_check_fasta(arabidopsis_proteome)
pombe_proteome, pombe_fasta = create_and_check_fasta(pombe_proteome)

# 2nd: create databases prior to performing reciprocal best BLAST. Specific databases types depend of the type of sequences in the file and would be created accordingly

make_db(pombe_fasta, pombe_proteome, dbname = "POMBE")
make_db(arabidopsis_fasta, arabidopsis_proteome, dbname = "ARABIDOPSIS")


# 3rd: Create factories to perform blast from own databases (just created)

ara_factory = Bio::Blast.local('blastp', './Databases/ARABIDOPSIS')
pombe_factory = Bio::Blast.local('blastp', './Databases/POMBE')


# 4th: Perform blast and parse the output to do best hits reciprocal blast to find putative othologues

putative_othologues_candidates = Hash.new   # hash to store the found candidates

# Each sequence of the S.pombe proteome would be blasted against the arabidopsis proteome
pombe_fasta.each_entry do |entry|
    # Perfom blast with each entry as the query sequence and retrieve the best hit if any
    first_blastp_best = blast_and_best_hit(entry, ara_factory)
    next if first_blastp_best.nil?  # if best hit not found, then pass to next sequence

    # Get original sequence from arabidopsis proteome from a given id to perfom the second blast
    first_blastp_best_id = first_blastp_best.target_def.split("|")[0].delete(' ')
    first_blastp_best_entry = search_sequence(first_blastp_best_id, arabidopsis_proteome)

    # Perform reciprocal blast and retrieve best hit
    second_blastp_best = blast_and_best_hit(first_blastp_best_entry, pombe_factory)    # Now the query is the best hit from Arabidopsis
    next if second_blastp_best.nil?
    
    second_blastp_best_id = second_blastp_best.target_def.split("|")[0].delete(' ') # retrieve id of the reciprocal best hit

    # If best hit from reciprocal blast matches the original S.pombe query (initial sequence or entry) then this and the first best hit are candidates to be putative orthologues
    if second_blastp_best_id == entry.entry_id.delete(' ') 
      puts "#{second_blastp_best_id} is an orthologue candidate to #{first_blastp_best_entry.entry_id}"
      putative_othologues_candidates[first_blastp_best_id] = second_blastp_best_id  # store pair of candidates found for printing the report
      puts putative_othologues_candidates.length
    else
        puts "-#{second_blastp_best_id}- is not equal to -#{entry.entry_id.delete(' ')}-"
    end
end

# Create a report with the retrieve orthologue candidates

output_report_file = './Orthologue_candidates_report.txt'
write_candiates_report(putative_othologues_candidates, output_report_file)