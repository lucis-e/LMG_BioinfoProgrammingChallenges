require 'rest-client'
require 'json'  # esto de momento no lo uso
require './NetworkMember'
require './InteractionNetwork'



# IntAct DB to retrieve protein - protein interaction
base_address = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/'

# Parameters
interactor_locus_name = 'At4g18960' # we should iterate in this, also we would first create a hash to change this
species = 'species:arabidopsis'  # filtering by species
tab25_format = 'tab25'  # retrieving in tab25 format

# Query url
query_adress = "#{base_address}search/interactor/#{interactor_locus_name}/?query=#{species}&format=#{tab25_format}"


# Get the response for our protein-protein information request 
response = RestClient::Request.execute(
    method: :get,
    url: query_adress,
    headers: {'Accept' => 'application/json'})  # use the RestClient::Request object's method "execute"
    

network = InteractionNetwork.new

network.read_intAct_response(response)


#-------------------------------------------------------------------------------------------------------
## PARA COMPROBAR QUE LA COSA FUNCIONA
para_comprobar_duplicados = []

network.network_members.keys.each do |instance|
  para_comprobar_duplicados << instance # append al array
  puts instance.other_ids
end

puts para_comprobar_duplicados.length
#puts 
#puts para_comprobar_duplicados.uniq! != nil   # uniq! modifies the array variable (removes duplicate elements) , returns nil if no changes were performed
# False =  no hay repetidos