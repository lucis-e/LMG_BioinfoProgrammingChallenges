require 'rest-client'
require 'json'  # esto de momento no lo uso
require './NetworkMember'
require './InteractionNetwork'


# IntAct DB to retrieve protein - protein interaction
# ----------------------------------------------------------
INTACT_BASE_ADDRESS = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/'

# Parameters for Intact DB
interactor_locus_name = 'AT4g18960' # we should iterate in this, also we would first create a hash to change this
SPECIES = 'species:arabidopsis'  # filtering by species
TAB25_FORMAT = 'tab25'  # retrieving in tab25 format

# Query url

#-----------------------------------------------------------------
# FUNCTION: Programatic access to REST APIs and get responses

def rest_api_request(adress)
  RestClient::Request.execute(
    method: :get,
    url: adress,
    headers: {'Accept' => 'application/json'})  # use the RestClient::Request object's method "execute"
end


#-----------------------------------------------------------------


#-----------------------------------------------------------------

INTERACTION_DEPTH = 2   # Constant, that is why it is in caps, depth = 2 so we look for direct interactors and direct interactos of our direct interactors

network = InteractionNetwork.new
network.get_interactions(interactor_locus_name, INTERACTION_DEPTH)




#-------------------------------------------------------------------------------------------------------
## PARA COMPROBAR QUE LA COSA FUNCIONA
para_comprobar_duplicados = []

network.network_members.keys.each do |instance|
  para_comprobar_duplicados << instance # append al array
  puts [instance.id_interactor, instance.gene_id]
end

puts para_comprobar_duplicados.length
#puts 
#puts para_comprobar_duplicados.uniq! != nil   # uniq! modifies the array variable (removes duplicate elements) , returns nil if no changes were performed
# False =  no hay repetidos

