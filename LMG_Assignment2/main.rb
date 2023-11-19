require 'rest-client'
require 'json'
require './Members'
require './Networks'

#------------------MAIN FUNCTIONS --------------------------------------------

def rest_api_request(adress)
  RestClient::Request.execute(
    method: :get,
    url: adress,
    headers: {'Accept' => 'application/json'})  # use the RestClient::Request object's method "execute"
end

def extract_xref(element) # syntax of TAB25 <XREF><VALUE>(<DESCRIPTION>)
  # match regex to get <VALUE>
  match_data = element.match(/:(\S+)(?:\(|$)/)  
  # get the actual value
  val = match_data[1] if match_data
  return(val)
end

#-----------------------------------------------------------------------------



Members.read_from_file('ArabidopsisSubNetwork_GeneList.txt')


# Parameters for this assignment
INTACT_BASE_ADDRESS = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/'
SPECIES = 'species:arabidopsis' 
TAB25_FORMAT = 'tab25'
DEPTH = 2


def recursive_search(member, network, depth)

  if depth < 1
    return "End of search"
  else
    protein_list.map{|element| element.find_interactors if member.direct_interactors.empty?}  # seach only when we have not already searched
    member.add_to_network
    if member.direct_interators.empty? || member.direct_interactors.is_a?(String)
      next  # go for next element in the array
    end 
    return recursive_search(member.direct_interactors, network, depth = depth + 1)
  end
  


end

Members.all_coexpresed_members.each do |gene|
  puts
  puts
  puts "ANALIZANDO NUEVO MIEMBRO con #{gene.gene_id} y #{gene.uniprot_id}"
  network = Networks.new()
  gene.set_network=(network)
  network.add_member(gene)
  recursive_search(gene, network, depth=1)
  puts "Este miembro tiene la red #{network} con miembros:"
  network.network_members.each do |_key, value|
    puts value.uniprot_id
  end
end