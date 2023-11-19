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


def recursive_search(member, network, depth)

  if depth == 1
    puts "Estamos en nivel #{depth}"
    member.find_interactors if member.direct_interactors.empty? 
    member.add_to_network
    if member.direct_interactors.empty? || member.direct_interactors.is_a?(String)
      puts "Interacciones para #{member.uniprot_id} no encontradas en IntAct, corte de la recursividad"
      return 1
    end
    return recursive_search(member.direct_interactors, network, depth = depth + 1)

  elsif depth == 2
    puts "Estamos en nivel #{depth}"
    list_of_interactors = []
    member.each do |one_member|
      one_member.find_interactors if one_member.direct_interactors.empty?
      one_member.add_to_network
      list_of_interactors += one_member.direct_interactors
    end

    return recursive_search(list_of_interactors, network, depth = depth + 1) 

  else
    puts "Hemos pasado el depth 2"
    return 1
  end

end

Members.all_coexpresed_members.each do |member|
  puts
  puts
  puts "ANALIZANDO NUEVO MIEMBRO con #{member.gene_id} y #{member.uniprot_id}"
  network = Networks.new()
  member.set_network=(network)
  network.add_member(member)
  recursive_search(member, network, depth=1)
  puts "Este miembro tiene la red #{network} con miembros:"
  network.network_members.each do |_key, value|
    puts value.uniprot_id
  end
end