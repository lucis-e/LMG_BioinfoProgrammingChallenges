require 'json'


class Members

    attr_accessor :uniprot_id, :times_searched, :direct_interactors

    @@coexpresed_members= []
    @@all_members = Hash.new

    def initialize( uniprot_id: "IDXXX" )
        @uniprot_id = uniprot_id
        @times_searched = 0
        @direct_interactors = []
        @@all_members[@uniprot_id] = self
    end

    def gene_id=(gene_name) #SET gene_id
        @gene_id = gene_name
    end

    def gene_id #GET gene_id
        @gene_id
    end

    def set_network=(network)
        @network = network
    end
    
    def get_network
        @network
    end

    def self.all_coexpresed_members
        @@coexpresed_members
    end

    def all_members
        @@all_members
    end
    
    def self.read_from_file(filename) # Leer archivo y crear Members con cada ATG...
        coexpressed_file = File.open(filename, 'r')
        coexpressed_file.readlines.each do |line|
            locus_name=line.chomp
            togo_address = "http://togows.dbcls.jp/entry/uniprot/#{locus_name}/accessions.json"
            togo_response = rest_api_request(togo_address)  # search in TOGO db the uniprot id for each locus name
            result = JSON.parse(togo_response.body)
            if result.is_a?(Array) && result.any?
                uniprot_id = result.first.first   
            else
                puts "No UniProt entry found for locus #{locus_name}. Please remove this entry from gene list"
                next
            end
            member = self.new(uniprot_id: uniprot_id)   # create new instance of this class for each gene of the list with uniprotid
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
            [0,1].each do |id|
                interactor = extract_xref(values[id])
                if !interactor.nil? && !interactor.include?(self.uniprot_id)    # check if it not empty or the query interactor
                    if @@all_members.include?(interactor) # Si ya existe
                        @direct_interactors << @@all_members[interactor]
                    else
                        interactor = self.class.new(uniprot_id: interactor) 
                        @direct_interactors << interactor
                    end # Si todavía no existe
                end
            end
        end
    end


    def add_to_network #Coge los @direct_interactors de la instancia Member y 1)Define su @network como el del objeto Member 2)Los añade al Networks del Member
        if @direct_interactors.empty? || @direct_interactors.is_a?(String)
            return 1
        end
        @direct_interactors.each do |interactor|
            unless self.get_network.network_members.key?(interactor.uniprot_id)
                interactor.set_network=(self.get_network)
                self.get_network.add_member(interactor)
            end
        end
    end

    def register_search
        @times_searched += 1
    end



    def eql?(other) 
        self.uniprot_id == other.uniprot_id if other.is_a?(Networks) # just to make sure we are comparing objects of the same class
    end

    def hash    
        @uniprot_id.hash
    end
end