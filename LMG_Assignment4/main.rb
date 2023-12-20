require 'bio'

# -----------------------CONSTANTS ----------------------- #

# Start and stop codons to verify CDS sequences
START_CODON = 'ATG'
STOP_CODON = ['TAA', 'TAG', 'TGA']
# Evalue threshold for blast report hits
EVALUE_THRESHOLD = 1e-6  

# ----------------------- MAIN FUNCTIONS ----------------------- #

# Checks fasta file type: nucleotide or protein
# @param fasta [String] fasta file to check
# @return [String] nucleotide or protein

def fasta_sequence(fasta)
    first_seq = fasta.first.seq
    biosequence = Bio::Sequence.new(first_seq).auto   # guess the sequence type of this file
    if biosequence.instance_of? Bio::Sequence::NA
        return 'nucl'
    elsif biosequence.instance_of? Bio::Sequence::AA
        return 'prot'
    end
end 

# Filter out those sequences that are not CDS
# @param nucleotide_fasta_file [String] the fasta file
# @param start_codon [String] the start codon
# @param stop_codons [Array] the stop codons
# @return [Array] all fasta entrys that are actually CDS

def only_cds?(nucleotide_fasta_file, start_codon, stop_codons)
    filtered_proteome = []
    nucleotide_fasta_file.each do |entry|
        next unless entry.seq.start_with?(start_codon)
        next unless stop_codons.any? {|stop_codon| entry.seq.end_with?(stop_codon)}
        filtered_proteome << entry
    end
    return filtered_proteome
end

# Translate fasta CDS file to fasta protein file
# @param cds_filtered [String] the cds file already filtered
# @param translated_filename [String] the name of the new file
# @return [Bio::Sequence] fasta file with proteins

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

# Creates Bio::FastaFormat file with proteins and checks only CDS sequences when nucleotide
# @param proteome_file [String] the file to create and check
# @return [proteome_file] fasta file name with proteins
# @return [proteome_fasta] fasta file with proteins

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


# Creates data base from fasta file
# @param fasta_file [Bio::FastaFormat] the file to create the data base
# @param filename [String] the file name
# @param dbname [String] the data base name
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

# Searchs best hit regarding e-value for a query blast
# @param entry [Bio::FastaFormat] the fasta entry as query
# @param factory [Bio::Blast] the blast factory to search
# @return [Bio::Blast] fthe best hit found for report

def blast_and_best_hit(entry, factory)
    query = ">#{entry.entry_id}\n#{entry.seq}"
    report = factory.query(query)
  
    filtered_hits = report.hits.select{|hit| !hit.evalue.nil? && hit.evalue <= EVALUE_THRESHOLD}
    sorted_hits = filtered_hits.sort_by {|hit| hit.evalue}

    best_hit = sorted_hits.first
    return best_hit
end

# Searchs for an entry of a file given a protein identifier
# @param identifier [String] the protein ID
# @param fasta_file [String] the fasta file where the protein is
# @return [Bio::FastaForma] the entry of with that ID

def search_sequence(identifier,fasta_file)
    fasta = Bio::FlatFile.open(Bio::FastaFormat, fasta_file)

    found_entry = fasta.find { |entry| entry.entry_id.match?(Regexp.escape(identifier))}
    if found_entry
        return found_entry
    else
        puts "COULD NOT FIND YOUR SEQUENCE"
        return
    end
end

# Writes report with putative orthologues genes
# @param cadidate_hash [Hash] the orthologues found
# @param output_report_file [String] the name of the file to report
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

# Check number of arguments for two files

if ARGV.length != 2
    abort "Incorrect number of files passed. Proteome files of the two species should be specified"
end

# Check order of arguments

if ARGV[0] == "TAIR10_cds.fa" && ARGV[1] == "proteome_pombe.fa"    
    arabidopsis_proteome = ARGV[0]
    pombe_proteome = ARGV[1]
elsif ARGV[1] == "TAIR10_cds.fa" && ARGV[0] == "proteome_pombe.fa"
    pombe_proteome = ARGV[0]
    arabidopsis_proteome = ARGV[1]
else
    abort "Incorrect files passed. Usage: ruby main.rb TAIR10_cds.fa proteome_pombe.fa"
end


# 1st: Create BioRuby objects from files and check them
# Translates nucleotide file to proteins if found

arabidopsis_proteome, arabidopsis_fasta = create_and_check_fasta(arabidopsis_proteome)
pombe_proteome, pombe_fasta = create_and_check_fasta(pombe_proteome)

# 2nd: Create both databases from files.

make_db(pombe_fasta, pombe_proteome, dbname = "POMBE")
make_db(arabidopsis_fasta, arabidopsis_proteome, dbname = "ARABIDOPSIS")


# 3rd: Create factories to perform blast search on both data bases

ara_factory = Bio::Blast.local('blastp', './Databases/ARABIDOPSIS')
pombe_factory = Bio::Blast.local('blastp', './Databases/POMBE')


# 4th: Perform best reciprocal hits blast search

# Creates a hash to store orthologues
putative_othologues_candidates = Hash.new

# Loops over one of the file entries
pombe_fasta.each_entry do |entry|
    # Retrieve first blastp best hit
    first_blastp_best = blast_and_best_hit(entry, ara_factory)
    next if first_blastp_best.nil?

    # Retrieve original sequence to perform second blast
    first_blastp_best_id = first_blastp_best.target_def.split("|")[0].delete(' ')
    first_blastp_best_entry = search_sequence(first_blastp_best_id, arabidopsis_proteome)

    # Retrieve second blastp best hit and its identifier
    second_blastp_best = blast_and_best_hit(first_blastp_best_entry, pombe_factory)    # perform second blastp with best hit from S.pombe
    next if second_blastp_best.nil?
    
    second_blastp_best_id = second_blastp_best.target_def.split("|")[0].delete(' ')

    # If second blastp best hit turns out to be the entry of the file parsed, then orthologues have
    # been found
    if second_blastp_best_id == entry.entry_id.delete(' ') 
      puts "#{second_blastp_best_id} is an orthologue candidate to #{first_blastp_best_entry.entry_id}"
      putative_othologues_candidates[first_blastp_best_id] = second_blastp_best_id
      puts putative_othologues_candidates.length
    else
        puts "-#{second_blastp_best_id}- is not equal to -#{entry.entry_id.delete(' ')}-"
    end
end

# Writes report with orthologue genes

output_report_file = './Orthologue_candidates_report.txt'
write_candiates_report(putative_othologues_candidates, output_report_file)