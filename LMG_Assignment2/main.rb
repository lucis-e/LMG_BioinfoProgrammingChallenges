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




Members.all_coexpresed_members.each do |gene|
  puts
  puts
  puts "ANALIZANDO NUEVO MIEMBRO con #{gene.gene_id} y #{gene.uniprot_id}"
  
  network = Networks.create_and_merge(gene, DEPTH)
  puts "Este miembro tiene la red #{network} con #{network.network_members.length} miembros"
  #network.network_members.each do |element|
  #  puts element.uniprot_id
  #end
end
  puts
  puts
  puts "FINAL REPORT"
  puts "--------------------------------------------------"
  puts "Total number of nets: #{Networks.get_number_of_nets}"
  

