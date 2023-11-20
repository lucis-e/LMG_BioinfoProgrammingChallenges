require 'rest-client'

INTACT_BASE_ADDRESS = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/'
SPECIES = 'species:arabidopsis' 
TAB25_FORMAT = 'tab25'


#prote = Q9M3B6
proteina = "P42737"

def find_interactors(prote, intact_address=INTACT_BASE_ADDRESS, species=SPECIES, formato=TAB25_FORMAT)
    intact_address = "#{intact_address}search/interactor/#{prote}/?query=#{species}&format=#{formato}"
    response = RestClient::Request.execute(
        method: :get,
        url: intact_address,
        headers: {'Accept' => 'application/json'})

    if response.empty? 
      puts  "Response Not Available in IntAct"
      return 1
    end
    response.body.each_line do |line|
        values = line.chomp.split("\t")

        #filtering by tecnoloogy
        interaction_detection_method = extract_xref(values[6])  # get confidence score 
        next if interaction_detection_method == "MI:0018"   # exclude two hybrid, not very trustworthy, a lot of false positives
        #puts interaction_detection_method

        # filering by type of interaction
        type_int = extract_xref(values[11])  # get confidence score
        next if type_int != "MI:0915"   # if it is not physical association then filter, reserved for purified molecules, not cell extracts or other complex sambples (100% no other participants are bridging the interaction)
        puts type_int
        #next if instact_miscore < 0.1 # filter, exclude those that have a confidence score lower than 0.3

    end
end

def extract_xref(element) # syntax of TAB25 <XREF><VALUE>(<DESCRIPTION>)
    # match regex to get <VALUE>
    match_data = element.match(/:(.*?)(?=\(|$)/)    
    # get the actual value
    val = match_data[1] if match_data
    return(val)
end


find_interactors(proteina)