class Networks

  attr_accessor :network_members

  @@number_of_networks = 0
  @@all_networks = Hash.new

  def initialize
    @network_members = {}
  end

  def add_member(net_member)
    @network_members[net_member.uniprot_id] = net_member
  end

  
  
end
