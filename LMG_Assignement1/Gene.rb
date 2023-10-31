# This code was created in collaboration with my colleague Miguel La Iglesia Mirones
# Lucía Muñoz Gil, MSc in Computational Biology, Bioinformatics Programming Challenges course, 2023
# Code for Assigment #1

# Code for defining the "Gene" Class Object

class Gene

    # Accesor to set/get the values of some of the properties from outside the Object
    attr_accessor :name, :mutant_phenotype

    # Initialize method defining some of the main properties and setting default values when an instance is created
    def initialize( 
        name: "XXXX",   
        mutant_phenotype: "NA") 
        @name = name    # Name of the gene, if no value is passed then name="XXXX"
        @mutant_phenotype = mutant_phenotype    # Mutant phenotype description, if no value is passed then is "NA"
    end

    # "Get" method for gene_ID property (to get the value from outside the Object)
    def gene_ID
        @gene_ID
    end

    # "Set" method for gene_ID property, the function is defined since we want to check for format
    def gene_ID=(string)
        unless string.match(/A[Tt]\d[Gg]\d\d\d\d\d/)    # Programm stops running if the Gene Identificator does not meet
                                                        # the required format. 
            warn "CODE STOP: Gene Identifier format is incorrect. It should have the format '/A[Tt]\d[Gg]\d\d\d\d\d/' "
            exit 1
        end
        @gene_ID = string   # If the Gene ID has the right format, then gene_ID property takes value
    end

    # "Get" method for linked_to property (to get the value of the property from outside the Object)
    def linked_to
        @linked_to
    end

    # "Set" method for linked_to property, the function is defined since only those "Gene" instances that are proved
    # to be genetically-linked to some other gene would have a value assigned to this property.
    def linked_to=(new_gen)
        if @linked_to   # Check if there is already a value assigned for this property (in case the gene is discovered to
                        # be linked to more than 1 gene).
            @linked_to += ",#{new_gen}" # If there is already a value then add the new gene name
        else
            @linked_to = new_gen    # If there is no value, then assign the gene name as value
        end
    end
end