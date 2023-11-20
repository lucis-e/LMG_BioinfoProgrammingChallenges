class Networks

  attr_accessor :network_members, :number_of_networks

  @@total_networks = 0
  @@all_networks = Hash.new

  def initialize
    @network_members = []
    @@total_networks += 1
    @@all_networks[self] = @network_members # esto habria que cambiarlo, no sé si es útil
  end

  def add_member(net_member)
    @network_members << net_member
  end

  def self.get_number_of_nets
    @@total_networks
  end

  def self.all_networks
    @@all_networks
  end


  def add_interactors_to_network(net_member)
    net_member.direct_interactors.each do |interactor|
        unless @network_members.include?(interactor)
            #interactor.set_network=(self)
            add_member(interactor)
        end
    end 
  end


  def recursive_search(found_proteins, depth)
    return if depth <= 0 # end searching for interactors when depth is 0
  
    list_of_interactors = []
    found_proteins.each do |protein|
      protein.find_interactors if protein.direct_interactors.empty? # search only if we have not already searched, if not use what is already in memory (info that we have already extracted, no need to search again!)
      next if protein.direct_interactors.empty? || protein.direct_interactors.is_a?(String)
  
      add_interactors_to_network(protein)
      list_of_interactors += protein.direct_interactors # use information just searched or that was already stored
    end
    
    recursive_search(list_of_interactors, depth - 1)
  end

  def self.create_and_merge(gene, depth)
    network = Networks.new  # create the new network
    #gene.set_network=(network)
    network.add_member(gene)
    network.recursive_search([gene], depth) # search and assign all the interactors found to this net

    existing_network == @@all_networks.values.find{ |network| network.network_members.keys.any?}

    return(network)
  end

end