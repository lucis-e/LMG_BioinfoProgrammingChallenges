
class Members

    attr_accessor :uniprot_id, :direct_interactors

    @@coexpresed_members= []
    @@all_members = Hash.new
    @@total_num_members = 0

    def initialize( uniprot_id: "IDXXX" )
        @uniprot_id = uniprot_id
        @direct_interactors = []
        @@all_members[@uniprot_id] = self
        @@total_num_members +=1
    end

    def gene_id=(gene_name) #SET gene_id
        @gene_id = gene_name
    end

    def gene_id #GET gene_id
        @gene_id
    end

    def self.all_coexpresed_members
        @@coexpresed_members
    end

    def self.all_members
        @@all_members
    end

    def self.number_of_members
        @@total_num_members
    end
    
    def self.read_from_file(filename) # Leer archivo y crear Members con cada ATG...
        coexpressed_file = File.open(filename, 'r')
        coexpressed_file.readlines.each do |line|
            locus_name=line.chomp

            if locus_name !~ /A[Tt]\d[Gg]\d\d\d\d\d/
                abort "Locus name #{locus_name} does not meet the correct format. Please define locus names as AT0g00000"
            end

            result = togo_search("uniprot", locus_name,"/accessions")
            if result.is_a?(Array) && result.any?
                uniprot_id = result.first.first   
            else
                puts "No UniProt entry found for locus #{locus_name}. Please remove this entry from gene list"
                next
            end
            member = AnnotatedMembers.new(uniprot_id: uniprot_id)   # create new instance of this class for each gene of the list with uniprotid
            member.gene_id=(locus_name) # and genename
            @@coexpresed_members << member
        end
    end

    

    def find_interactors(intact_address=INTACT_BASE_ADDRESS, species=SPECIES, formato=TAB25_FORMAT)
        intact_address = "#{intact_address}search/interactor/#{@uniprot_id}/?query=#{species}&format=#{formato}"
        response = rest_api_request(intact_address)
        if response.empty? 
          @direct_interactors = "Response Not Available in IntAct"
          return 1
        end
        response.body.each_line do |line|
            values = line.chomp.split("\t")

            # Filtering by confidence score
            instact_miscore = extract_xref(values[14]).to_f   # get confidence score 
            next if instact_miscore < 0.4 # filter, exclude those that have a confidence score lower than 0.35
            
            # Filtering by quality of/trust in tecnology
            interaction_detection_method = extract_xref(values[6])  # get detection method
            next if interaction_detection_method == "MI:0018"   # exclude two hybrid, not very trustworthy, a lot of false positives

            # Filtering by type of interaction
            type_int = extract_xref(values[11])  # get type of interaction
            next if type_int != "MI:0915" # if it is not physical association then filter, reserved for purified molecules, not cell extracts or other complex sambples (100% no other participants are bridging the interaction)
                    
            [0,1].each do |id|
                interactor = extract_xref(values[id])
                if !interactor.nil? && !interactor.include?(@uniprot_id) && interactor.match(/[OPQ][0-9][A-Z0-9]{3}[0-9]$/)   # check if it not empty or the query interactor
                    if @@all_members.key?(interactor) # Si ya existe
                        @direct_interactors << @@all_members[interactor]
                    else
                        interactor = AnnotatedMembers.new(uniprot_id: interactor) 
                        @direct_interactors << interactor
                    end 
                end
            end
        end
    end

    def eql?(other) 
        self.uniprot_id == other.uniprot_id if other.is_a?(Networks) 
    end

    def hash    
        @uniprot_id.hash
    end

end