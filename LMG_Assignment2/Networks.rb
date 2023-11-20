class Networks

  attr_accessor :network_members, :number_of_networks

  @@total_networks = 0
  @@all_networks = []

  def initialize
    @network_members = []
    @@total_networks += 1
    @@all_networks << self # más optimo un array
  end

  def self.all_networks
    return @@all_networks
  end


  def add_member(net_member)
    @network_members << net_member
  end

  def self.get_number_of_nets
    @@total_networks
  end

  def only_one_member?
    @network_members.length < 2
  end

  def self.reduce_networks
    @@all_networks = @@all_networks.select { |network| !network.only_one_member? }
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


  def merge_with(other_network)

    common_members = @network_members & other_network.network_members # check if there are any common members (this is actually already checked, maybe not necessary)

    if common_members.any?
      @network_members |= other_network.network_members # |= simulates "union", now the @network_members would have all unique members from both nets
      @@all_networks.delete(other_network)  # deleting the old net, now we have a new one with its info and new info
      @@total_networks -= 1 # decrease in 1 the total number of nets, we have just merge 2 into 1
    end
  end



  def create_and_merge
                                                                        # Lu hace lo de ir absorbiendo redes igual q yo, por eso el select te selecciona aquellso que cumplen la condiicon de abajoo
    nets_with_common_members = @@all_networks.select do |existing_net|  # iterate through all existing nets y selecciona aquellos que cumplan esa condición
      existing_net.network_members.any? { |member| self.network_members.include?(member)}  # check if there is any net with common members with the just created net
    end # El any? mira si alguno de los elementos del @network_members cumple la condición de {"se incluye en la red del objeto que llama a la función"}
      # it is possible that the new net has common members with more than 1 net, we should merge them all in that case

    if nets_with_common_members.size > 1
      nets_with_common_members.each do |common_net|  
        self.merge_with(common_net) unless common_net == self   # if would match to itself, we need to skip because if not we will be deleting the own net we are creating
      end
    end

    return(nets_with_common_members)
  end

end