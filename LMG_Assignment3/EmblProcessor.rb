class EmblProcessor
    
    attr_accessor :entries

    # the emblprocessor is the embl flat file so entries are instances of the class entry
    def initialize(filename)
        @entries = scan_repetitive_features(filename)
    end

    def scan_repetitive_features(filename)
        embl_file = Bio::FlatFile.auto(filename)    # crea el objeto embl flat file
        entries = []

        embl_file.each_entry do |entry| # como antes recorre las entrys
            next unless entry.accession # Lack of accesion is suspicious
            start_entry, end_entry = get_coordinates(entry.definition)  #  Extract sequence chromosome coordinates from entry.definition    
            embl_entry = EmblEntry.new(entry, start_entry, end_entry)   # New EmblEntry instance with: entry and chromosomal coordinates
            embl_entry.annotate_source_gene
            entries << embl_entry   # add each entry instance to the array with all the entrys
        end
        
        return entries
    end  

    def load_to_gff(coordinates_to_use)
        gff = Bio::GFF::GFF3.new
        count = 1
    
        @entries.each do |embl_entry|
          embl_entry.process_cttctt_repeats(coordinates_to_use)
          
          feature_exists = embl_entry.entry_sequence.features.any? { |feature| feature.feature == 'cttctt_repeat' }
          
          if !feature_exists 
            puts "#{embl_entry.gene_locus} has no exons with the CTTCTT repeat" if coordinates_to_use != "chromosomal coordinates"
            next
          end

          embl_entry.entry_sequence.features.each do |feature|

            next unless feature.feature == 'cttctt_repeat'
    
            qual = feature.assoc
            attributes = [{ 'ID' => "repeat_region_#{count}", 'Name' => "cttctt_repeat_#{embl_entry.gene_locus}" }]
    
            seqid = embl_entry.seq_id
    
            gff.records << Bio::GFF::GFF3::Record.new(
              seqid,            # seqID
              'programmatically',   # source
              qual['SO_Name'],  # feature type
              qual['start'],    # start
              qual['end'],      # end
              nil,              # score
              qual['strand'],   # strand
              nil,              # phase
              attributes[0]     # attributes
            )
            count += 1 
          end
        end
    
        return gff
      end


    
end