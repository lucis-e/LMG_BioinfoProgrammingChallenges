class AnnotatedMembers < Members

    attr_accessor :go_IDs_terms, :kegg_gene, :kegg_ID_pathway
    
    def initialize(**args) #(go_IDs:, go_terms:, kegg_gene:, kegg_ID_pathway:)
        super
        @go_IDs_terms = {}
        @kegg_gene = ""
        @kegg_ID_pathway = {}
        self.annotate_GO
        self.annotate_kegg
    end

    def annotate_kegg
        return if @kegg_gene.empty?

        result = togo_search("kegg-genes", @kegg_gene)
        if result[0].key?("pathways") && !result[0]["pathways"].nil? && !result[0]["pathways"].empty? 
            @kegg_ID_pathway = result[0]["pathways"]
        end
        #puts "EL hash es: #{@kegg_ID_pathway}"
    end


    def annotate_GO
        result = togo_search("ebi-uniprot", self.uniprot_id, "/dr")
        if result[0].key?("GO") && !result[0]["GO"].nil? && !result[0]["GO"].empty?
            list_of_GOs = result[0]["GO"]
            list_of_GOs.each do |go|
                biological_process = go[1].match(/^[P]:/) #molecular functions (F), biological processes (P), and cellular components (C).
                if biological_process
                    id = go[0]
                    term = biological_process.post_match
                    @go_IDs_terms[id] = term
                else
                    next
                end
            end
        end

        if result[0].key?("KEGG") && !result[0]["KEGG"][0][0].nil? && !result[0]["KEGG"][0][0].empty?
            @kegg_gene = result[0]["KEGG"][0][0]
        end
        
    end

end