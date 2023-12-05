#------------------------ PUBLIC INSTANCE METHODS -----------------------------------------------------------------------

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
    return start_coordinate, end_coordinate
end


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


def scan_repetitive_features(filename, coordinates_to_use)

    embl_file = Bio::FlatFile.auto(filename)

    search_positive = Bio::Sequence::NA.new("cttctt")
    search_complementary = Bio::Sequence::NA.new("aagaag")

    regex_positive = Regexp.new(search_positive.to_re)
    regex_complementary = Regexp.new(search_complementary.to_re)

    i = 0

    source = "" # esto es necesario definirlo así pues si no tiene problemas
    entries = []
    k=0
   

    # Loop each entry in EMBL file
    embl_file.each_entry do |entry|

        #break if i > 20
        #i += 1

        next unless entry.accession # Lack of accesion is suspicious

        # Extract sequence chromosome coordinates from entry.definition
        start_entry, end_entry = get_coordinates(entry.definition)

        # Anotate sequence chromosome coordinates
        entry_sequence = entry.to_biosequence

        f1 = Bio::Feature.new('chromosomal_coordinates', entry.accession) # 'feature type' , 'position'
        f1.append(Bio::Feature::Qualifier.new('start', start_entry))
        f1.append(Bio::Feature::Qualifier.new('end', end_entry))
        entry_sequence.features << f1

        # LOOP FEATURES
        entry.features.each do |feature|
            featuretype = feature.feature  
            position = feature.position 
            qual = feature.assoc

            # Anotate the source for feature
            if featuretype == 'source'
                source = qual['db_xref']
            end

            # To annotate Arabidopsis thaliana gene as a new feature. Tests if it was not annotated before and if the format of coordinates are correct.
            if featuretype == 'gene' && feature.position.match(/^(complement\()?(\d+\.\.\d+)(\))?$/)

                # Test if that feature already exists
                feature_exists = entry_sequence.features.any? { |feature| feature.feature == 'gene_from_list'}
                next if feature_exists

                # New feature with AT gene name ("AT5G48300")
                f1 = Bio::Feature.new('gene_from_list', qual['gene']) # 'feature type' , 'position'
                entry_sequence.features << f1
            end
            
            # Tests if feature is "exon" and if the coordinates do not correspond to another fragment
            unless featuretype == 'exon' && position.match?(/^(complement\()?(\d+\.\.\d+)(\))?$/)
                next
            end

            # Get exon ID
            exon_id = qual["note"]

            # Get exon coordinates
            start_exon, end_exon = get_coordinates(position)
            
            # Get exon sequence from its coordinates, and test if that exon is inside the entry sequence to continue or not
            exon_sequence = entry_sequence.subseq(start_exon, end_exon)
            next if exon_sequence.nil? || exon_sequence.empty?

            # Decides wether to use positive or complement regular expression
            if position.include?('complement')
                regex = regex_complementary
                strand = '-'
                coordinates_format = "complement(%s)" # %s indicates the substituient, in this case, coordinates range
            else
                regex = regex_positive
                strand = '+'
                coordinates_format = "%s"
            end

            # Find motifs in exon sequence and retrieves the positions (BUT IN EXON SEQUENCE). It works for several matches/motifs in the exon.
            matches_motifs = find_motifs(exon_sequence, regex)
            
            next if matches_motifs.empty? || matches_motifs.nil?
            
            # For each matched motif in exon sequence
            matches_motifs.each do |match|

                # First and last position of motif, BUT IN EXON SEQUENCE
                start_match = match[:start]
                end_match = match[:end]

                # Conversion from exon positions to entry sequence positions
                start_motif = start_exon + start_match # +1 por ccs del match, -1 por la suma de exon y match
                end_motif = start_exon + end_match - 1 # -1 por la suma de exon y match
                #puts "Motif #{entry_sequence.subseq(start_motif, end_motif)} found at position #{start_motif} to #{end_motif} in #{exon_id}"
            
                # If chromosomal coordinates, it recalculates those motif coordinates in the chromosomal coordinates
                if coordinates_to_use == "chromosomal coordinates"
                    chr_coordinates = entry_sequence.features.find { |feature| feature.feature == 'chromosomal_coordinates' }
                    qual_chr = chr_coordinates.assoc
                    start_chr = qual_chr['start']
                    end_chr = qual_chr['end']

                    length_motif = end_motif - start_motif
                    start_motif = start_motif + start_chr - 1
                    end_motif = start_motif + length_motif - 1 # PROBAR CON EL EXON_SEQUENCE?    
                    
                end
                

                coordinates = coordinates_format % "#{start_motif}..#{end_motif}"

                # Tests if that feature has already been annotated, as several exons have same ctttcttt positions
                feature_exists = entry_sequence.features.any? { |feature| feature.feature == 'cttctt_repeat' && feature.position == coordinates }
                next if feature_exists

                k = k + 1

                # Anotates new feature
                f1 = Bio::Feature.new('cttctt_repeat', coordinates ) # 'feature type' , 'position'
                f1.append(Bio::Feature::Qualifier.new('start', start_motif))
                f1.append(Bio::Feature::Qualifier.new('end', end_motif))
                f1.append(Bio::Feature::Qualifier.new('strand', strand))
                #f1.append(Bio::Feature::Qualifier.new('score', '.')) # In this case, there is no score value
                #f1.append(Bio::Feature::Qualifier.new('phase', '.')) # This is only set when feature type is "CDS", not "exon"
                f1.append(Bio::Feature::Qualifier.new('source', source))
                f1.append(Bio::Feature::Qualifier.new('SO_Name', 'repeat_region'))
                entry_sequence.features << f1


            end

        end

        # Keeping all entry objects in a list for later
        entries << entry_sequence

    end
    return entries
end


def load_to_gff(entries, coordinates_to_use)
    gff = Bio::GFF::GFF3.new
    h=0
    count = 1 # To annotate attributes

    entries.each do |entry|

        # To test if that entry has a repetitive motif annotated, this value True/False will be used later
        feature_exists = entry.features.any? { |feature| feature.feature == 'cttctt_repeat' }

        entry.features.each do |feature|

            featuretype = feature.feature
            qual = feature.assoc
            
            # For those without repetitive motif, puts it gene locus ("AT5G48300"). We avoid repetition when chromosomal coordinates
            if featuretype == 'gene_from_list' && !feature_exists && coordinates_to_use != "chromosomal coordinates"
                # ACTIVATE THIS PIECE OF CODE
                puts "#{feature.position} has no exons with the CTTCTT repeat"
                break
            end

            next unless featuretype == 'cttctt_repeat'
            h = h + 1
            
            #puts "Chr#{entry.entry_id} #{qual['source'].class} #{qual['SO_Name'].class} #{qual['start'].class} #{qual['end'].class} #{qual['strand'].class} ID=repeat_region_#{count};Name=CTTCTT_motif"
            
            attributes = [{"ID" => "repeat_region_#{count}", "Name" => "CTTCTT_motif"}]
            
            # If chromosomal coordinates, seqid field gives chromosomal coordinates of the entry
            if coordinates_to_use == "chromosomal coordinates"
                chr_coordinates = entry.features.find { |feature| feature.feature == 'chromosomal_coordinates' }
                seqid= chr_coordinates.position
            else
                seqid = "Chr#{entry.entry_id}"
            end

            gff.records << Bio::GFF::GFF3::Record.new(
                seqid,      # seqID
                qual['source'],     # source
                qual['SO_Name'],    # feature type
                qual['start'],            # start
                qual['end'],          # end
                nil,          # score
                qual['strand'],          # strand
                nil,          # phase
                attributes[0] # attributes
            )
            count += 1
        end
    end
    return gff
end

def write_gff(gff, coordinates_to_use)
    if coordinates_to_use == "chromosomal coordinates"
        file = 'AT_repeats_chromosomal.gff'
    else
        file = 'AT_repeats_sequence.gff'
    end 

    File.open(file, 'w') do |file|
        file.puts gff.to_s
    end
end