# == NetworkMember
# This blbalblabla
#
# == Summary
  

class NetworkMember

    attr_accessor :id_interactor, :other_ids

# Get/Set the new 
# @!attribute [rw]
# @return [String] The namethi

    def initialize(
        id_interactor:"IDXXX",
        other_ids:"IDXXX"
    )
        @id_interactor = id_interactor  # ojo hay dos members por interaccion, es decir linea de la reques
        @other_ids = other_ids

    end

    # Set method of gene id, not every member would have gene id, only interactors in the gene lists
    def gene_id=(gene_name)
        @gene_id = gene_name
    end

    # Get method for gene_id, this would be useful for the report I imagine
    def gene_id
        @gene_id
    end
    
    # This would be made so intances with the same @id_interactor would be considered as equal (independent from object_id)
    # This is useful when working with hashes, we ensure that objects with the same content are considered equañ
    
    # key? uses eql? mehtod to compare? and that is why we have to override so it considers two objects equals if having the same value for id_interactor attribute
    def eql?(other) 
        self.id_interactor == other.id_interactor if other.is_a?(NetworkMember) # just to make sure we are comparing objects of the same class
    end

    def hash    # this code generates a hash code based on the attribute values
        # it is important for the correct fucntioning of hash-based collections (like Ruby's Hash), since we are storing our netmembers in a hash, we do this so when we look for duplicates, we do it by id
        @id_interactor.hash
    end
end