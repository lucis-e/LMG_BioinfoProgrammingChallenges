module LuMike_objects

    # == Members
    #
    # Nodes of protein-protein interaction networks. Members have protein ID (Uniprot ID), direct interactors (if any) and gene locus name (any)
    #
    # == Summary
    # 
    # Class used for repressenting protein members of a protein-protein interaction network. Members are proteins that have an unique identifier, also some also can have the name of the locus encoding for the protein.
    # If direct interactors of a member are found, they are recorded as an attribute of this member, and the interactors become themselves Members of the name network
    #

    class Members

        # Get/Set method for the Uniprot ID of the protein
        # @!attribute [rw] uniprot_id
        # @return [String] Protein Uniprot ID
        attr_accessor :uniprot_id

        # Get/Set mtehod for the direct interactors found to be binding the protein
        # @!attribute [rw] direct_interactors
        # @return [array<Members>] List of direct interactos (also instances of the class Members)
        attr_accessor :direct_interactors
        
        # @atrr [Array<Members>] Members which are provided as input genes (input file) and supossed to be coexpressed
        @@coexpresed_members= []

        # @atrr [Hash] All instances of the Members class, accessed by key Uniprot ID, even if not part of an specific network
        @@all_members = Hash.new

        # @atrr [Integer] Number of members (nodes), even if not part of a network
        @@total_num_members = 0

        # Create a new instance of the class Members 
        # @param uniprot_id [String] Uniprot ID of the protein
        # @return [Members] Just created member
        def initialize( uniprot_id: "IDXXX" )
            # @atrr [String] Protein Uniprot ID
            @uniprot_id = uniprot_id 
            # @atrr [Array<Members>] List of direct interactos 
            @direct_interactors = []                
            @@all_members[@uniprot_id] = self       
            @@total_num_members +=1                 
        end

        # Set method for locus name of the gene encoding for the protein member
        # @param gene_name [String] Gene locus name 
        # @return [void]
        def gene_id=(gene_name)
            @gene_id = gene_name    
        end
        
        # Get the gene_id of the member
        # @!method no_params_method
        # @return [String] the gene ID of the member
        def gene_id
            @gene_id
        end
        
        # Get method for retrieving all members inputed as coexpressed genes
        # @!method no_params_method
        # @return [Array<Members>] List of "coexpressed" genes
        def self.all_coexpresed_members
            @@coexpresed_members
        end

        # Get method for retrieving all Members instances created
        # @!method no_params_method
        # @return [Array<Members>] List of all existing nodes/members
        def self.all_members
            @@all_members
        end

        # Get method for retrieving the total number of members
        # @!method no_params_method
        # @return [Integer] Total number of members
        def self.number_of_members
            @@total_num_members
        end

        # Process the input file of "coexpressed" genes, create AnnotatedMember instances for each gene if meting the required format
        # Prior, make request to TOGO for retrieving the Uniprot ID of the protein each specific gene is coding for
        # The AnnotatedMember is annotated with GO Tems and KEGG Pathways
        # Add member to list of coexpressed genes
        # @param filename [String] Filename of inputed file with the supposedly coexpressed genes 
        # @return [void]
        def self.read_from_file(filename)
            coexpressed_file = File.open(filename, 'r')
            coexpressed_file.readlines.each do |line|
                locus_name=line.chomp

                if locus_name !~ /A[Tt]\d[Gg]\d\d\d\d\d/ # Check for correct format of gene locus name
                    abort "Locus name #{locus_name} does not meet the correct format. Please define locus names as AT0g00000"
                end

                result = togo_search("uniprot", locus_name,"/accessions") # Request TOGO database to retrieve the Uniprot ID for the encoded protein
                if result.is_a?(Array) && result.any?
                    uniprot_id = result.first.first 
                else    # Fitering: we are only taking into account those genes registered in Uniprotkb
                    puts "No UniProt entry found for locus #{locus_name}. Please remove this entry from gene list"
                    next
                end
                member = LuMike_objects::AnnotatedMembers.new(uniprot_id: uniprot_id) 
                member.gene_id=(locus_name)
                @@coexpresed_members << member
            end
        end

        
        
        # Search for interaction information in Intact database given a member. 
        # Make REST API request for accesing Intact database and retrieve protein-protein interaction info of the query protein that meet the quality requirements
        # Process response (two proteins per interaction) and create and store annotated members if not already recorded in previous searches
        # @param intact_address [String] Base Intact address for constructing the final query URL
        # @param species [String] Species proteins belong to (filtering)
        # @param formato [String] Desired output format
        # @return [void]
        def find_interactors(intact_address=INTACT_BASE_ADDRESS, species=SPECIES, formato=TAB25_FORMAT)
            intact_address = "#{intact_address}search/interactor/#{@uniprot_id}/?query=#{species}&format=#{formato}" # Access through PSICQUIC REST API
            response = rest_api_request(intact_address)
            if response.empty? 
            @direct_interactors = "Response Not Available in IntAct"
            return 1
            end
            response.body.each_line do |line|
                values = line.chomp.split("\t")

                # FILTER: confidence score. We assume a 0.5 is sufficient for a trustworthy interaction
                instact_miscore = extract_xref(values[14]).to_f
                next if instact_miscore < 0.5
                
                # FILTER: innteraction detection tecnology. Exclude those interactions detected by two-hybrid technology (high rate of false positives)
                interaction_detection_method = extract_xref(values[6])
                next if interaction_detection_method == "MI:0018"

                # FILTER: type of interaction should be direct interaction (reserved for purified molecules, not complex molecules). 100% sure not other participants are enabling of the interaction
                type_int = extract_xref(values[11])
                next if type_int != "MI:0915"
                
                [0,1].each do |id| 
                    interactor = extract_xref(values[id])   # A Member would only be created if it meets the requirements and it has not been already created (or is the query member)
                    if !interactor.nil? && !interactor.include?(@uniprot_id) && interactor.match(/[OPQ][0-9][A-Z0-9]{3}[0-9]$/) # Doublecheck: we'd only select those interactors with a standard uniprot id
                        if @@all_members.key?(interactor) 
                            @direct_interactors << @@all_members[interactor]
                        else
                            interactor = LuMike_objects::AnnotatedMembers.new(uniprot_id: interactor) 
                            @direct_interactors << interactor   # record as interactor of the query member
                        end 
                    end
                end
            end
        end
    end
end