module LuMike_objects

    # == AnnotatedMembers
    #
    # Nodes of protein-protein interaction networks (Members) once they have annotated with the pathways these members are part of (KEGG Pathways) and its related GO Terms
    #
    # == Summary
    # 
    # Class used for repressenting members of protein-protein interaction networks considering their molecular characteristics in the cell
    #

    class AnnotatedMembers < LuMike_objects::Members

        # Get/Set method for the GO ID and Term associated to the member
        # @!attribute [rw]
        # @return [Hash] Keys are the GO IDs and values are the GO Terms
        attr_accessor :go_IDs_terms

        # Get/Set method for the identificator of KEGG associated with this member's locus name
        # @!attribute [rw]
        # @return [String] KEGG id for member's locus name
        attr_accessor :kegg_gene

        # Get/Set method for the KEGG ID and pathways the member is part of
        # @!attribute [rw]
        # @return [Hash] Keys are the KEGG IDs and values are the pathway descriptions
        attr_accessor :kegg_ID_pathway

        # Create a new instance of the class AnnotatedMembers
        # @return [AnnotatedMembers] Just created annotated member (inherits attributes and methods from Members)
        def initialize(**args)
            super
            # @atrr [Hash] GO ID and Terms associated to the member/protein
            @go_IDs_terms = {}
            # @atrr [String] KEGG identificator for the member's locus name
            @kegg_gene = ""
            # @#atrr [Hash] KEGG ID and description of the molecular pathway the member is part of 
            @kegg_ID_pathway = {}
            self.annotate_GO
            self.annotate_kegg
        end

        # Given the KEGG id for the locus name, make a TOGO REST API request to search in KEGG database for the ID and the description of the molecular pathways the protein is part of 
        # @!method no_params_method
        # @return [void]
        def annotate_kegg
            return if @kegg_gene.empty? # Check for presence of KEGG id, required to search for the molecular pathway details

            result = togo_search("kegg-genes", @kegg_gene)
            if result[0].key?("pathways") && !result[0]["pathways"].nil? && !result[0]["pathways"].empty? 
                @kegg_ID_pathway = result[0]["pathways"]
            end
        end
        
        # Given the Uniprot ID of the member, make a TOGO REST API request to search in GO database for the GO IDs and Go Terms associated to this member/protein
        # Addtionaly, it retrieves also information related to the kegg id for the locus name of the member, if any (for KEGG annotations)
        # Filtering: record only those GO Terms related to biological processes
        # @!method no_params_method
        # @return [void]
        def annotate_GO
            result = togo_search("ebi-uniprot", self.uniprot_id, "/dr")
            if result[0].key?("GO") && !result[0]["GO"].nil? && !result[0]["GO"].empty? 
                list_of_GOs = result[0]["GO"]   # Only GO terms would be recorded
                list_of_GOs.each do |go|
                    biological_process = go[1].match(/^[P]:/) # Filering: Biological processes (P)
                    if biological_process
                        id = go[0]
                        term = biological_process.post_match    # What is after the GO ID is the GO Term
                        @go_IDs_terms[id] = term
                    else
                        next
                    end
                end
            end

            if result[0].key?("KEGG") && !result[0]["KEGG"][0][0].nil? && !result[0]["KEGG"][0][0].empty?   # Retrieve the KEGG ID associated to the locus name if any
                @kegg_gene = result[0]["KEGG"][0][0]
            end  
        end
    end
end