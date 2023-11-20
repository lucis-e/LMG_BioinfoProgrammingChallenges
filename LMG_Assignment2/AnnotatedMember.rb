class AnnotatedMember < Member
    #COGER LO DE LA CARPETA DE "IMPORTANTE"
    attr_accessor :kegg_id, :pathway_name, :go_id, :go_term
    
    def initialize
        super

    end

    def annotate_kegg

    end

    def annotate_GO

    end

end