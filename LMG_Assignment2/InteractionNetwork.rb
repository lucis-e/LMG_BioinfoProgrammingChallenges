# == InteractionNetwork
#
#
#
# == Summary
#

class InteractionNetwork

  attr_accessor :network_members   # no tengo claro poner accesor, o solo un get method

  @@number_of_networks = 0    # we would like a count of how many nets
  @@all_network_object = Hash.new   # registry of all instances of interactionnetwork object 


  def initialize
    @network_members = {}  # Lets create an emty dictionary for indexing the network members (yo diria de meter el objeto aqui, o sea la instancia, ero no tengo claro esto de crear el diccionario vacío, por lo de mergear)
  end


  # Add member to interactionnetwork, we should put some id
  def add_member(net_member)
    @network_members[net_member] = true   # asi decimos que el objeto está en el array, se comprueba de forma más sencilla
  end


  def extract_xref(element) # syntax es <XREF><VALUE>(<DESCRIPTION>)
    # match regex to get <VALUE>
    # match regex to get the value 
    match_data = element.match(/:(\S+)(?:\(|$)/)  #HAY QUE EXPLICAR ESTO BN!! Que menudo lio, no termino de entenderlo pero hay que verlo
    # get the actual value
    val = match_data[1] if match_data
    return(val)
  end

  # Togo REST API to get the Uniprot from locus name: IntAct querys we would like to do them with UniProt (protein interactions, not interactions with gene)
  def get_UniprotID_from_TOGO(locus_name)
    togo_uri = "http://togows.dbcls.jp/entry/uniprot/#{locus_name}/accessions.json"

    togo_response = rest_api_request(togo_uri)  # Make request to Togo REST API
    result = JSON.parse(togo_response.body)  # parse json into array
    
    if result.is_a?(Array) && result.any? # check if we have retrieved any response
      uniprot_id = result.first.first   # get the first uniprot id from list, there might be more than one
      return uniprot_id
    else
      puts "No UniProt entry found for locus #{locus_name}. Please remove this entry from gene list"
      return nil
    end

  end


  # ESTO O LO DIVIDIMOS EN FUNCIONES PEQUEÑAS O CREAMOS OTRA CLASE INTERACTIONS.rb y LO HACEMOS PORQUE ESTO ES UN LIO Y ADEMÑAS HAY QUE METERLE MÁS COSAS  esta es la que tiene que ser recursia
  
  def get_interactions(interactor_locus_name, depth, intact_address = INTACT_BASE_ADDRESS, species=SPECIES, format_intact=TAB25_FORMAT)

    if interactor_locus_name.match(/AT\dg\d\d\d\d\d/)
      @gene_id = interactor_locus_name  # save gene id, we would set this as attribute in member class
      interactor_uniprot = get_UniprotID_from_TOGO(interactor_locus_name)
      # va a haber que crear una clase para los genes igual, con el id de prote y el id de gen
      # o creamos la instncia y la metemos en hashhh
      # hostia que buena yo creo que esa esta bn, porque siiiii la creamos y la metemos aqui en el hash directamete
      # luego en el unless de abajo, no la va a meter y yo estaba pensando poner la funcion recursiva depues del unles
      # o sea, llamar a la misma función despues del unless, o sea ya solo con los nuevo sque meta en el hash, para no 
      # estar buscando varias veces los mismos, entonces si el inicial ya lo metemos, no va a pasar el corte en el unless
      #porque ya va a estar y no va a volver a buscarlo (entoces quito lo de los otros ids y yas ta, si no croe que los usemos)
    else
      interactor_uniprot = interactor_locus_name  # depth 2, looking for interactors of the interactors
    end

    query_adress = "#{intact_address}search/interactor/#{interactor_uniprot}/?query=#{species}&format=#{format_intact}"

    # Get the response for our protein-protein information request (IntAct)
    intAct_response = rest_api_request(query_adress)

    intAct_response.body.each_line do |line|
      values = line.chomp.split("\t") # split each line into array (tab25 with 15 columns/fields, everyone means something)

      # Get detection_method for this interaction (filtering?)
      interaction_detection_method = extract_xref(values[6])  # por si queremos filtrar por esto (MI_0045 es experimental interaction pero esto no está, hay que crear un hash con todos o ver que onda)

      # Type of interaction: esto es lo de physical association (MI:0915)
      type_of_int = extract_xref(values[11])  # esto es más facil de filtrar peroo no sé si menos riguroso

      # intact_miscore: por si queremos filtrar
      intact_miscore = values[13].to_f  # es un numero decimal, convert to float

      # una de las dos siempre va a ser la propia proteina (gen) que estamos buscando
      # tendremos que crear su clase una vez y las demas veces pues entiendo que no
      # aquí tendríamos que crear una clase:
      [0, 1].each do |id| # we would have to create 2 objects since the query protein can be A or B

        # Create the member
        member = NetworkMember.new(
          id_interactor:extract_xref(values[id]),  # if 0 it would be interactor A, if 1 it would be interactor, extract the value out of the xref syntax
          other_ids:extract_xref(values[id+2])  # get other id, values[2] if 0 for member A, values[3] if 1 for member b        
        )
        
        # Check if we already recorded this member as part of this network (we would also have to check for other networks) (lo hacemos más adelante, esto es el mergin)
        
        unless @network_members.key?(member) # if we don't use the unless there woyld not be duplicates also but because the objects would keep overwrittng, lets skip overwriting the same object again and again
          add_member(member) 
        end

      end

    end

  end

end




