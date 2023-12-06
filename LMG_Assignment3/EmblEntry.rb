
# == EmblEntry
#
# Object representing each entry of the .embl file (gene's sequence informatiion)
#
# == Summary
#
# Class used to represent information about the sequences of a given gene or list of genes, employing objects and methods of BioRuby gem
#

class EmblEntry

    # Get/Set method for the genomic sequence of a given gene
    # @!attribute [rw] entry_sequence
    # return [Bio::Sequence] the genomic sequence 
    attr_accessor :entry_sequence

    # Get/Set method for the chromosomal coordinates of a given gene
    # @!attribute [rw] chromosomal_coordinates
    # return [Array<Integer>] the chromosomal coordinates as a two-elemnt array
    attr_accessor :chromosomal_coordinates

    # Get/Set method for the source of the annotation for a given gene
    # @!attribute [rw] source
    # return [String] source of annotation
    attr_accessor :source

    # Get/Set method for the accession ID of a gene/entry
    # @!attribute [rw] accession
    # return [String] accession ID
    attr_accessor :accession

    # Get/Set method for the sequence ID of a gene/entry
    # @!attribute [rw] seq_id
    # return [Array<Integer>] the sequence ID
    attr_accessor :seq_id

    # Get/Set method for the gene locus name of a given gene
    # @!attribute [rw] gene_locus
    # return [Array<Integer>] the gene locus name
    attr_accessor :gene_locus

    # It creates a new instance of EmblEntry
    #
    # @param entry[Bio::EMBL] the BioRuby EMBL entry of a given gene
    # @param start_entry [Integer] the start chromosomal coordinate of the gene
    # @param end_entry [Integer] the en chromosomal coordinate of the gene
    # return [EmblEntry] Just created EMBL entry
    def initialize(entry, start_entry, end_entry)
        @accession = entry.accession  # get accesion ID from entry
        @seq_id = entry.entry_id
        @entry_sequence = entry.to_biosequence  # get nucleotide sequence of the gene
        @chromosomal_coordinates = [start_entry, end_entry] # create an array with the chromosomic coordinates of the gene
        @source = ""
        @gene_locus =""
    end

    # It annotates the gene locus of the gene as a feature of the Bio::Sequence (genomic sequence) object
    # @!method no_params_method
    # return [void]
    def annotate_source_gene
        @entry_sequence.features.each do |feature|
          featuretype = feature.feature
          position = feature.position
          qual = feature.assoc
          
          # Annotate source
          if featuretype == 'source'
            @source = qual['db_xref']
          end
          
          # Annotate gene locus for the sequence
          if featuretype == 'gene' && position.match(/^(complement\()?(\d+\.\.\d+)(\))?$/) && @gene_locus.empty?
            @gene_locus = qual['gene']
          end

        end
        return
    end


    # Adds a found CTTCTT motif as feature of the EMBL Bio::Sequence object
    #
    # @param coordinates [String] the start..end coordinates of the motif (in the gene's sequence or chromosome)
    # @param start_motif [Integer] the start coordinate of the motif (gene's sequence or chromosomal)
    # @param end_motif [Integer] the end coordinate of the motif (gene's sequence or chromosomal)
    # @param strand [String] genomic strand were motif can be found (+ or -)
    # return [void]
    def add_cttctt_repeat(coordinates, start_motif, end_motif, strand)
        return if @entry_sequence.features.any? { |f| f.feature == 'cttctt_repeat' && f.position == coordinates }   # check if this motif has already been annotated
    
        f1 = Bio::Feature.new('cttctt_repeat', coordinates) # New feature of the Embl biosequence object --> repetitive motifs "CTTCTT"
        f1.append(Bio::Feature::Qualifier.new('start', start_motif))
        f1.append(Bio::Feature::Qualifier.new('end', end_motif))
        f1.append(Bio::Feature::Qualifier.new('strand', strand))
        f1.append(Bio::Feature::Qualifier.new('SO_Name', 'repeat_region'))
        @entry_sequence.features << f1  # record as feature
    end

    # Converts a motif gene's sequence coordenates to chromosomal coordenates
    #
    # @param start_motif [Integer] the start coordinate of the motif (gene's sequence or chromosomal)
    # @param end_motif [Integer] the end coordinate of the motif (gene's sequence or chromosomal)
    # return [Array<Integer>] the chromosomal coordinates of the CTTCTT motif
    def convert_to_chromosomal_coordinates(start_motif, end_motif)
        start_chr = @chromosomal_coordinates[0]
        length_motif = end_motif - start_motif   # length of the motif
        start_motif = start_motif + start_chr - 1 # - 1 so 1 + 1 is still 1
        end_motif = start_motif + length_motif   # end coordinate of the motif in the chromosome
        return [start_motif, end_motif]
    end

    # Finds CTTCTT repeats in exons in either the positive strand or the complement. It retrieves coordinates dependent of the reference system (chromosome or sequence)
    #
    # @param coordinates_to_use [String] reference system to calculate motif's coordinates
    # @param regex_positive [Regexp] regular expression for searching for the CTTCTT pattern in the positive strand
    # @param regex_complementary [Regeexp] regular expresion for searching for the CTTCTT pattern in the complementary strand 
    # return [void]
    def process_cttctt_repeats(coordinates_to_use, regex_positive=REGEX_POSITIVE, regex_complementary = REGEX_COMPLEMENT)

        @entry_sequence.features.each do |feature|
          next unless feature.feature == 'exon' && feature.position.match?(/^(complement\()?(\d+\.\.\d+)(\))?$/)    # only search in exons
    
          start_exon, end_exon = get_coordinates(feature.position)  # retrive exonic coordenates
          exon_sequence = @entry_sequence.subseq(start_exon, end_exon)  # get exon subsequence from exon coordinates
    
          next if exon_sequence.nil? || exon_sequence.empty?    # control: is empty?
    
          if feature.position.include?('complement')    # get directionality of exon (positive or complement strand)
            regex = regex_complementary
            strand = '-'
            coordinates_format = "complement(%s)" # specify the format in which the coordenates are given 
          else
            regex = regex_positive
            strand = '+'
            coordinates_format = "%s"
          end
        
          
          matches_motifs = find_motifs(exon_sequence, regex)  # find motifs and calculate its coordenates if any (exon is the reference system). Works for multiple motifs
          next if matches_motifs.empty? || matches_motifs.nil?
        
          # For each matched motif in exon sequence
          matches_motifs.each do |match|
            start_match = match[:start] # get exonic coordenates of motif
            end_match = match[:end]
            # Conversion from exon positions to entry sequence positions
            start_motif = start_exon + start_match  # +1 for object index in ruby, -1 because of the sum of coordenates
            end_motif = start_exon + end_match - 1  # -1 because of the sum of coodenates
    
             # If chromosomal coordinates, it recalculates those motif coordinates in the chromosomal coordinates
            if coordinates_to_use == 'chromosomal coordinates'
                start_motif, end_motif = convert_to_chromosomal_coordinates(start_motif, end_motif)
            end

            coordinates = coordinates_format % "#{start_motif}..#{end_motif}"
            add_cttctt_repeat(coordinates, start_motif, end_motif, strand, @source) # add motif as feature of the Bio::Sequence object
            
          end
        end
      end

end