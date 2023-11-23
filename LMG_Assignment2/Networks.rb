module LuMike_objects

  # == Networks
  #
  # Protein-protein interaction networks, consisting of members (nodes) and the GO and KEGG Pathways annotations related to its members
  #
  # == Summary
  # 
  # Networks of molecular interactors (proteins) that are known to bind eachother, becoming part of regulatory networks. Each node has some interactors which the member directly binds to 
  #


  class Networks

    # Get/Set mtehod for the list of members forming oart of interaction network
    # @!attribute [rw]
    # @return [Array<Members>] All Members instances part of this network
    attr_accessor :network_members

    # Get/Set mtheod for the total number of networks
    # @!attribute [rw]
    # @return [Integer] Total number of networks
    attr_accessor :number_of_networks

    # Get/Set method for the network GO annotations (related to its members)
    # @!attribute [rw]
    # @return [Hash] Union of the GO IDs and Terms associated to the members of the network
    attr_accessor :network_GOs

    # Get/Set method for the network KEGG annotations (related to its members)
    # @!attribute [rw]
    # @return [Hash] Union of the molecular pathways' IDs and descriptions the members of the network are involved in
    attr_accessor :network_KEGGs

    # @atrr [Integer] Total number of networks
    @@total_networks = 0
    # @atrr [Integer] Total number of "coexpressed" genes (genes described in input file) that do not belong to any network
    @@nodes_without_net = 0
    # @atrr [Array<Networks>] All instances of the class Networks ever created
    @@all_networks = []
    
    # Create a new instance of Networks
    # @!method no_params_method
    # @return [Networks] Just created network
    def initialize
      # @atrr [array<Members>] All Members instances part of this network
      @network_members = []     
      # @atrr [Integer] Total number of networks
      @@total_networks += 1
      # @atrr [Hash] Union of the GO IDs and Terms associated to the members of the network
      @network_GOs = {}
      # @atrr [KEGG] Union of the molecular pathways' IDs and descriptions the members of the network are involved in
      @network_KEGGs = {}
      @@all_networks << self
    end

    # Get method for retrieving all Networks instances created
    # @!method no_params_method
    # @return [array<Networks>] List of all networks
    def self.all_networks
      return @@all_networks
    end
    
    # Get method for retrieving the total number of networks
    # @!method no_params_method
    # @return [Integer] Total number of networks
    def self.get_number_of_nets
      @@total_networks
    end

    # Get method for retrieving number of initial gene members not part of any network
    # @!method no_params_method
    # @return [array<Members>] Total number of "coexpressed" members not part of any network
    def self.nodes_without_net
      @@nodes_without_net
    end

    # Add a new member as part of the network
    # @param net_member [<Members>] Member to be added
    # @return [void]
    def add_member(net_member)
      @network_members << net_member
    end

    # Check for networks with less than 1 member (we do not consider isolated nodes as networks)
    # @!method no_params_method
    # @return [void]
    def only_one_member?
      @network_members.length < 2
    end

    # After all nets have been created and merged, delete network instances with 1 member from all_networks and reduce total number of networks
    # Addtionally, the number of excluded networks would be the number of isolated nodes not associated with any net (nodes_without_net)
    # @!method no_params_method
    # @return [void]
    def self.reduce_networks
      before_reduction = @@all_networks.length
      @@all_networks = @@all_networks.select { |network| !network.only_one_member? }

      @@nodes_without_net = before_reduction - @@all_networks.length  # number of "coexpressed" genes not part of any network
      @@total_networks -= @@nodes_without_net # reduce total number of nets
    end

    # Annotate the network. Records the union of the GO and KEGG annotations for the members that form the network
    # @!method no_params_method
    # @return [void]
    def annotate_network

      @network_members.each do |member| 

        # Annotate network with Go Terms associated to its members (only if they had not been already recorded)
        if !member.go_IDs_terms.nil? && !member.go_IDs_terms.empty?
          member.go_IDs_terms.each do |key, value|
            @network_GOs[key] = value unless @network_GOs.key?(key)
          end
        end

        # Annotate network with KEGG Pathways associated to its members (only if they had not been already recorded)
        if !member.kegg_ID_pathway.nil? && !member.kegg_ID_pathway.empty?
          member.kegg_ID_pathway.each do |key, value|
            @network_KEGGs[key] = value unless @network_KEGGs.key?(key) 
          end
        end

      end
    end


    # Add the direct interactors of a member of the network as part of this same network, only if they are not already part of this net
    # @param net_member [<Members>] Member already part of the net 
    # @return [void]
    def add_interactors_to_network(net_member)
      net_member.direct_interactors.each do |interactor|
        unless @network_members.include?(interactor) 
              add_member(interactor)
        end
      end 
    end


    # Recursive search for interactors in a network
    # @param found_proteins [Array<Members] the list of members to search for its interactors in the network
    # @return [void]
    def recursive_search(found_proteins, depth)
      return if depth <= 0
    
      list_of_interactors = []
      found_proteins.each do |protein|
        protein.find_interactors if protein.direct_interactors.empty?
        next if protein.direct_interactors.empty? || protein.direct_interactors.is_a?(String)
    
        add_interactors_to_network(protein)
        list_of_interactors += protein.direct_interactors # Interactors are stored for calling recursivity later
      end
      
      recursive_search(list_of_interactors, depth - 1)
    end

    # Merges two networks into one. It assigns the union of the members to the self network and deletes the other network, reducir total_networs by one
    # @param other_network [<Networks>] Networks object with common members to self
    # @return [void]
    def merge_with(other_network)

      common_members = @network_members & other_network.network_members # Tests if there are members in common
      if common_members.any?
        @network_members |= other_network.network_members # merge members into self network
        @@all_networks.delete(other_network) 
        @@total_networks -= 1 
      end
    end

    # Search for common members of the network with any other already created network. If common members are found then employs the merge_with method for merging all networks with common members to the network (becoming just one big net)
    # @!method no_params_method
    # @return [void]
    def merge_with_common

      nets_with_common_members = @@all_networks.select do |existing_net|  # check if there are networks with common members with self
        existing_net.network_members.any? { |member| self.network_members.include?(member)}
      end

      if nets_with_common_members.size > 1  # since the net would match to itself, merge only if there are more than 1 matches
        nets_with_common_members.each do |common_net|
          self.merge_with(common_net) unless common_net == self # merge networks two by two but not merge self with itself
        end
      end
    end
  end
end