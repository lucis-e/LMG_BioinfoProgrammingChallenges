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

# Ask user to input from comand line if wants to translate the CDS file for the proteome of some species
def ask_for_translation(file)
    loop do
        puts "Do you want to translate #{file}? (yes/no):"
        response = $stdin.gets.chomp.downcase

        case response
        when "no", "n"
            return false
        when "yes", "y"
            return true
        else 
            puts "Invalid response. Please enter yes (y) or no (n)"
        end
    end
end


# Translate CDS to protein: in case that user wants to perform reciprocal blastp

def translate_cds_to_protein(cds_file, translated_filename)
    proteome_file = File.open(translated_filename, 'w')
    cds_file.each_entry do |entry|
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
        if ask_for_translation(proteome_file)
            proteome_file = "TAIR10_translated.fa"
            proteome_fasta = translate_cds_to_protein(proteome_fasta, proteome_file)
        end
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
# Prompt the user to specify from command line which type of search to do
# translate (y/n): if yes then a reciprocal blastp is performed, if not then tblastn + blastx

arabidopsis_proteome, arabidopsis_fasta = create_and_check_fasta(arabidopsis_proteome)
pombe_proteome, pombe_fasta = create_and_check_fasta(pombe_proteome)


puts only_cds?(arabidopsis_fasta, START_CODON, STOP_CODON)

# 2nd:reate both databases prior to performing reciprocal best BLAST. Specific databases types depend of the type of sequences in the file and would be created accordingly

make_db(pombe_fasta, pombe_proteome, dbname = "POMBE")
make_db(arabidopsis_fasta, arabidopsis_proteome, dbname = "ARABIDOPSIS")


# 3rd: Create factories to perform blast depending on what kind of blast the user wants to perform

ara_factory = Bio::Blast.local('tblastn', './Databases/ARABIDOPSIS')
pombe_factory = Bio::Blast.local('blastx', './Databases/POMBE')


# 4th: Perform blast and parse the output to do best hits reciprocal blast to find putative othologues

putative_othologues_candidates = Hash.new

arabidopsis_fasta.each_entry do |entry|
    query = ">#{entry.entry_id}\n#{entry.seq}"
    report = pombe_factory.query(query)
  
    tblastn_best_hit = nil
  
    report.each_hit do |hit|
      if !hit.evalue.nil? && hit.evalue <= EVALUE_THRESHOLD
        tblastn_best_hit = hit if tblastn_best_hit.nil? || hit.evalue < tblastn_best_hit.evalue
        puts "#{hit.hit_id} : evalue #{hit.evalue}\t#{hit.target_id} "
      end
    end
    puts tblastn_best_hit
  
    next if tblastn_best_hit.nil?
  
    query = ">#{tblastn_best_hit.hit_id}\n#{tblastn_best_hit.target_seq}"
    puts query
    report = ara_factory.query(query)
  
    blastx_best_hit = nil
  
    report.each_hit do |hit|
      if !hit.evalue.nil? && hit.evalue <= EVALUE_THRESHOLD
  
        blastx_best_hit = hit if blastx_best_hit.nil? || hit.evalue < blastx_best_hit.evalue
        #puts "#{hit.hit_id} : evalue #{hit.evalue}\t#{hit.target_id} "
      end
    end
    #puts blastx_best_hit
  
    next if blastx_best_hit.nil?
  
    if blastx_best_hit.hit_id == entry.entry_id # Si el blastx de la proteina->CDS da el ID del entry-CDS de TAIR
      puts "#{entry.entry_id} is an orthologue candidate to #{tblastn_best_hit.hit_id}"
      putative_othologues_candidates[entry.entry_id] = tblastn_best_hit.hit_id
    end
end


def write_candiates_report()

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
    file.puts
    file.puts "---------------------------------------------------------------------------------------"
    file.puts
    end
end 