# == EmblProcessor
#
# Object representing a processor for .embl files to create a final GFF3 file
#
# == Summary
# 
# Class used for processing an EMBL file, creating objects for each of its entries (information about the sequence of a given gene) as EMBL objects and finally creating a GFF3 file with the annotated CTTCTT motifs found in exonic regions
#

class EmblProcessor
    
    # Get/Set method for all EmblEntry objects created (all entries of the .embl file)
    # @!attribute [rw]
    # @return [String] all EmblEntry instances
    attr_accessor :entries

    
    # It creates a new instance of EmblProcessor
    #
    # @param filename [String] the name of the .embl file
    # @return [EmblProcessor] Just created EMBL Processor
    def initialize(filename)
        @entries = scan_entries(filename)
    end

    # Create EmblEntry objects out of the entries of the EMBL file 
    # @param filename [String] the name of the .embl file
    # @return [entries] all EmblEntry instances created out of the EMBL file
    def scan_entries(filename)
        embl_file = Bio::FlatFile.auto(filename)  # .embl file as a Bio::FlatFile object 
        entries = []

        embl_file.each_entry do |entry|
            next unless entry.accession # If not accession, then next 
            start_entry, end_entry = get_coordinates(entry.definition)  #  Get chromosomal coordinates of each gene
            embl_entry = EmblEntry.new(entry, start_entry, end_entry) # create a new EmblEntry instance
            embl_entry.annotate_source_gene # Add as features the source and the gene locus name
            entries << embl_entry   
        end
        
        return entries
    end  

    # Creates a Bio::GFF:GFF3 object with all motifs annotated found in each gene exons to be printed out
    # @param coordinates_to_use [String] reference system to calculate the coordinates of the motif with respect to
    # @return [Bio::GFF::GFF3] the Bio::GFF::GFF3 object
    def load_to_gff(coordinates_to_use)
        
        report_count = 0 # ID of genes with no CTTCTT motif found in its exons
        entries_without_regions = []

        gff = Bio::GFF::GFF3.new  # new Bio::GFF::GFF3 object
        feature_count = 1 # ID of feature to be annotated in GFF3
    
        @entries.each do |embl_entry|

          embl_entry.process_cttctt_repeats(coordinates_to_use) # found all motifs of all exons and add it as feature of the entry /gene
          
          feature_exists = embl_entry.entry_sequence.features.any? { |feature| feature.feature == 'cttctt_repeat' } # check: has this entry/gene have been annotated
          
          if !feature_exists  # if not then this gene does not have CTTCTT patterns in its exons
            entries_without_regions << embl_entry
            report_count += 1
            next
          end

          embl_entry.entry_sequence.features.each do |feature|

            next unless feature.feature == 'cttctt_repeat'  # look only for CTTCTT motif features
    
            qual = feature.assoc
            attributes = [{ 'ID' => "repeat_region_#{feature_count}", 'Name' => "cttctt_repeat_#{embl_entry.gene_locus}" }] # atrributes of the GFF3 file
    
            seqid = embl_entry.seq_id
            
            # New record for GFF file
            gff.records << Bio::GFF::GFF3::Record.new(
              seqid,            # seqID
              'programatically',# source
              qual['SO_Name'],  # feature type
              qual['start'],    # start
              qual['end'],      # end
              nil,              # score
              qual['strand'],   # strand
              nil,              # phase
              attributes[0]     # attributes
            )
            feature_count += 1
          end
        end

        # REPORT: number and locus name of genes without CTTCTT pattern in its exons
        report = File.open("gff_report.txt", "w")
        report.puts "REPORT OF GENES WITHOUT \'CTTCTT\' REPETITIVE REGIONS AFTER CREATING GFF FILE"
        report.puts "---------------------------------------------------------------------------"
        entries_without_regions.each do |entry|
          report.puts "#{entry.gene_locus}"
        end
        report.puts
        report.puts "TOTAL GENES: #{report_count}"
        report.close

        return gff  # final Bio::GFF::GFF3 object
      end


    
end