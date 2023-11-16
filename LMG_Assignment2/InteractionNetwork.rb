# == InteractionNetwork
#
#
#
# == Summary
#

class InteractionNetwork

  attr_accessor :network_members   # no tengo claro poner accesor, o solo un get method

  @@number_of_networks = 0

  # Get/Set the patient's name
  # @!attribute [rw]
  # @return [String] The name
  
  def initialize
    @network_members = {}  # Lets create an emty dictionary for indexing the network members (yo diria de meter el objeto aqui, o sea la instancia, ero no tengo claro esto de crear el diccionario vacío, por lo de mergear)
  end

  # Add member to Iteraction member
  # @!method?
  # Adds NetworkMember class object to hash/array (no me decido)
  def add_member(net_member)
    @network_members[net_member] = true   # asi decimos que el objeto está en el array, se comprueba de forma más sencilla
  end

  def extract_xref(element) # syntax es <XREF><VALUE>(<DESCRIPTION>)
    # match regex to get <VALUE>
    
    # match regex to get the value 
    match_data = element.match(/:(\S+)(?:\(|$)/)

    # get the actual value
    val = match_data[1] if match_data

    return(val)
  end

  # Creo que esto puede servir para mergear instancias, habría uqe crear una nueva y mergear las dos
  def merge(another_interactionnetwork)
    @net_members.merge!(another_interactionnetwork.network_members)
  end

  
  def read_intAct_response(intAct_response)

    intAct_response.body.each_line do |line|
      values = line.chomp.split("\t") # split each line into array (tab25 with 15 columns/fields, everyone means something)

      # Get detection_method for this interaction (filtering?)
      interaction_detection_method = extract_xref(values[6])  # por si queremos filtrar por esto (MI_0045 es experimental interaction pero esto no está, hay que crear un hash con todos o ver que onda)

      # Type of interaction: esto es lo de physical association (MI:0915)
      type_of_int = extract_xref(values[11])  # esto es más facil de filtrar peroo no sé si menos riguroso

      # intact_miscore
      intact_miscore = values[13].to_f  # es un numero decimal, convert to float

      # una de las dos siempre va a ser la propia proteina (gen) que estamos buscando
      # tendremos que crear su clase una vez y las demas veces pues entiendo que no
      # aquí tendríamos que crear una clase:
      [0, 1].each do |id| # we would have to create 2 objects since the query protein can be A or B

        # Create the member
        member = NetworkMember.new(
          id_interactor:extract_xref(values[id]),  # if 0 it would be interactor A, if 1 it would be interactor, extract the value out of the xref syntax
          other_ids:extract_xref(values[id+2]),  # get other id, values[2] if 0 for member A, values[3] if 1 for member b
          member_alias:extract_xref(values[id+4]) # same logic for extracting both alias depending on what interactor are we creating the instance for
        )

        # Check if we already recorded this member as part of this network (we would also have to check for other networks) (lo hacemos más adelante, esto es el mergin)
        
        add_member(member) unless @network_members.key?(member) # record as member of the net unless it is already

      end

    end

  end

end




