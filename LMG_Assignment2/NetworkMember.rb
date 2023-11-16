# == NetworkMember
# This blbalblabla
#
# == Summary
  

class NetworkMember

    attr_accessor :id_interactor, :other_ids, :member_alias

# Get/Set the new 
# @!attribute [rw]
# @return [String] The namethi

    def initialize(
        id_interactor:"IDXXX",
        other_ids:"IDXXX",
        member_alias:'aaa'
    )
        @id_interactor = id_interactor  # ojo hay dos members por interaccion, es decir linea de la reques
        @other_ids = other_ids
        @alias = member_alias
    end

    def ==(other) # this function would evaluate if the new instance is equal to any other in the network_members BY ID_INTERACTOR 
        @id_interactor == other.id_interactor   
    end

    def hash    # this code generates a hash code based on the attribute values
        # it is important for the correct fucntioning of hash-based collections (like Ruby's Hash), since we are storing our netmembers in a hash, we do this so when we look for duplicates, we do it by id
        @id_interactor.hash
    end
    
    def eql?(other) # eql method is calling the == funtion, by these we are ensuring hash uses the == method when checking for equality in the hash
        self==other # compare self to elements in hash
    end
end