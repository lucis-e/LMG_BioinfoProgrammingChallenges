# REVISAR NCBI LISTA IDs
# exon_sequence = entry_sequence.subseq(start_exon, end_exon): si el exon no cae en seq es que no está?
# COMO EVITAR OVERLAP DE MISMAS FEATURES
# REVISAR QUE TODAS LAS LAS FEATURES ANOTADAS TENGAN UN FORMATO ESTÁNDAR EN LA PAGINA ESA
# EL ID EN ATTRIBUTES ES COMO QUIERAS PERO IGUAL PARA TODOS Y CON UNA ENUMERACIÓN NO?ç
# Ver SI EL SOURCE ES ESO O EL "annotated by Araport11", yo creo que eso
# DOMINIOS CTTCTT QUE SE REFIEREN A LA MISMA SECUENCIA, PERO QUE ESTÁN EN 2 EXONES DISTINTOS, CÓMO CONSIDERARLO?
# UNA SOLA LINEA DEL GFF Y LUEGO EN ATRIBUTES LE METES TODOS LOS EXONES QUE LO TIENEN?
#Añadir una linea de metadata??
# hacer un filtro por "note" del exon que incluya el gen de la lista ? realmente no es necesario pues los de ccs de la secuencia son los q se miraran y por
#tanto son del gen de la lista.
#CONFIRMAR SOURCE Y ATTRIBUTES y ARREGLAR LO DE "source = "" AL ppio"



#----------------- LOAD EXTERNAL LIBRARIES AND MODULES-------------------------------------------------------------------

require 'bio'
require 'net/http'
require './my_functions'

#------------------------------ INPUT FILES NAMES AS ARGUMENTS FROM  COMMAND LINE ---------------------------------------

# Check for one common error: not specifying the correct number of files needed for the program to run
if ARGV.length != 1
    abort "Incorrect number of files passed. Two files names must be specified: input list of genes and name for final_report"
end
  
# Check for second common error: incorrect usage, files in incorrect order or wrong name passed.
#if ARGV[0] == "ArabidopsisSubNetwork_GeneList.txt"
if ARGV[0] == 'soloungen.txt'
        gene_file = ARGV[0]
else 
    abort "Incorrect order of files passed. Usage: ruby main.rb ArabidopsisSubNetwork_GeneList.txt"
end



#-------------------------------------------MAIN CODE ------------------------------------------------------------------


# BEFORE EVERYTHING, SEARCHING EMBL FILE (Task 1)

        # 1:  Using BioRuby, examine the sequences of the ~167 Arabidopsis genes from the last assignment by retrieving them from whatever database you wish #
        # ArabidopsisSubNetwork file -> each gene locus (lines) -> search for embl file -> add to 'AT_sequences.embl'

puts "Processing #{gene_file} file, this might take a while..."

        # Create string with gene IDs from file
gene_ids = read_from_file(gene_file)

        # Only one search with list of IDs
response_body = ncbi_fetch(database = 'ensemblgenomesgene',file_format = 'embl', id = gene_ids)

output_file = File.open('prueba.embl', 'w')
output_file.write(response_body)
output_file.close

# FIRST PART WITH SEQUENCE COORDINATES (Tasks 2, 3, 4a, 4b)

coordinates_to_use = "sequences coordinates"
entries = scan_repetitive_features('prueba.embl', coordinates_to_use)
gff = load_to_gff(entries, coordinates_to_use)
write_gff(gff, coordinates_to_use)

# SECOND PART WITH CHROMOSOMAL COORDINATES (Tasks 5, 6)

coordinates_to_use = "chromosomal coordinates"
entries = scan_repetitive_features('prueba.embl', coordinates_to_use)
gff = load_to_gff(entries, coordinates_to_use)
write_gff(gff, coordinates_to_use)