require 'bio'

# cosas que se le pueden añadir
# Parametros para hacer el blast
# comprobar qeu son cds
# que hacemso depsues?
# opcion traducir o no desde linea de comandos o como parametro del script
# que decida que tipo de blast en función del tipo de secuencia

# -----------------------CONSTANTS ----------------------------------------------

START_CODON = 'ATG'
STOP_CODON = ['TAA', 'TAG', 'TGA']
EVALUE_THRESHOLD = 1e-6  # Threshold e-value as proposed in # meter ref que no me deja copiar y pegar en la maquina virtual


# -------------- MAIN FUNCTIONS -------------------------------------------------

# Check content of the fasta files of the proteomes, the sequences are genomic sequences or protein sequences
def fasta_sequence(fasta)
    first_seq = fasta.first.seq
    biosequence = Bio::Sequence.new(first_seq).auto   # guess the sequence type of this file
    if biosequence.instance_of? Bio::Sequence::NA
        return 'nucl'
    elsif biosequence.instance_of? Bio::Sequence::AA
        return 'prot'
    end
end 

# Filter out those sequences that do not start with a start codon and end with an end codo

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


# CREATE AND CHECK FASTA
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
def make_db(fasta_file, filename, dbname)
    sequence_type = fasta_sequence(fasta_file)    # get the sequence type to define -dbtype
    result = system("makeblastdb -in #{filename} -dbtype '#{sequence_type}' -out ./Databases/#{dbname}")
    if $?.success?
        puts "BLAST database #{dbname} created successfully.\n"
    else
        puts "Error creating BLAST database. Output:\n#{result}\n"
    end
end

# Perform blast and get best hit based on e-value
def blast_and_best_hit(entry, factory)
    query = ">#{entry.entry_id}\n#{entry.seq}"
    report = factory.query(query)
  
    filtered_hits = report.hits.select{|hit| !hit.evalue.nil? && hit.evalue <= EVALUE_THRESHOLD}
    sorted_hits = filtered_hits.sort_by {|hit| hit.evalue}

    best_hit = sorted_hits.first
    return best_hit
end

# Search for the original sequence in the proteome to perform second blast (and not use the result of the alignment)
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

# Check input: introducing the file names of the proteomes of the species to discover putativo orthologues among 
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


# 1st: create BioRuby objects and check for the sequence type of the fasta files. 
# Translate the proteome file with cds nucleotide sequences into proteins to perform reciprocal blastp

arabidopsis_proteome, arabidopsis_fasta = create_and_check_fasta(arabidopsis_proteome)
pombe_proteome, pombe_fasta = create_and_check_fasta(pombe_proteome)

# 2nd:create both databases prior to performing reciprocal best BLAST. Specific databases types depend of the type of sequences in the file and would be created accordingly

make_db(pombe_fasta, pombe_proteome, dbname = "POMBE")
make_db(arabidopsis_fasta, arabidopsis_proteome, dbname = "ARABIDOPSIS")


# 3rd: Create factories to perform blast depending on what kind of blast the user wants to perform

ara_factory = Bio::Blast.local('blastp', './Databases/ARABIDOPSIS')
pombe_factory = Bio::Blast.local('blastp', './Databases/POMBE')


# 4th: Perform blast and parse the output to do best hits reciprocal blast to find putative othologues

putative_othologues_candidates = Hash.new

pombe_fasta.each_entry do |entry|
    first_blastp_best = blast_and_best_hit(entry, ara_factory)    # perform first blastp with protein A from Arabidopsis
    next if first_blastp_best.nil?

    # Retrieve original sequence to perform second blast
    first_blastp_best_id = first_blastp_best.target_def.split("|")[0].delete(' ')
    first_blastp_best_entry = search_sequence(first_blastp_best_id, arabidopsis_proteome)

    second_blastp_best = blast_and_best_hit(first_blastp_best_entry, pombe_factory)    # perform second blastp with best hit from S.pombe
    next if second_blastp_best.nil?
    
    second_blastp_best_id = second_blastp_best.target_def.split("|")[0].delete(' ')

    if second_blastp_best_id == entry.entry_id.delete(' ') 
      puts "#{second_blastp_best_id} is an orthologue candidate to #{first_blastp_best_entry.entry_id}"
      putative_othologues_candidates[first_blastp_best_id] = second_blastp_best_id
      puts putative_othologues_candidates.length
    else
        puts "-#{second_blastp_best_id}- is not equal to -#{entry.entry_id.delete(' ')}-"
    end
end

output_report_file = './Orthologue_candidates_report.txt'

write_candiates_report(putative_othologues_candidates, output_report_file)