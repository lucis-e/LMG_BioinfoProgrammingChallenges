require 'rest-client'
require 'json'
require './Members'
require './Networks'
require './AnnotatedMembers'

#------------------MAIN FUNCTIONS --------------------------------------------

def rest_api_request(address)
  RestClient::Request.execute(
    method: :get,
    url: address,
    headers: {'Accept' => 'application/json'})  # use the RestClient::Request object's method "execute"
end

def extract_xref(element) # syntax of TAB25 <XREF><VALUE>(<DESCRIPTION>)
  # match regex to get <VALUE>
  match_data = element.match(/:(\S+)(?:\(|$)/)  # la he cambiado
  # get the actual value
  val = match_data[1].gsub(/\A"/, '').gsub(/"\z/, '') if match_data  # get the value if any and if it has quotes then remove them
  return(val)
end

def togo_search(database,query,field="")
  togo_address = "http://togows.dbcls.jp/entry/#{database}/#{query}#{field}.json"
  togo_response = rest_api_request(togo_address)
  result = JSON.parse(togo_response.body)
  return result
end

#---- INPUT FILE NAMES AS ARGUMENTS FROM COMMAND LINE -------------------------------------------

# First warning
if ARGV.length != 2 # Check for one common human error: not specifying all 4 file names needed for the assignement
  abort "Incorrect number of files passed. Two files names must be specified: input list of genes and name for final_report"
end

# Second warning: Check for another common human error: not specifying the files in the correct order
if ARGV[0] == "ArabidopsisSubNetwork_GeneList.txt" && ARGV[1] == "Final_report.txt"
  input_gene_list = ARGV[0]
  output_report_file = ARGV[1]
else 
  abort "Incorrect order of files passed. Usage: ruby main.rb ArabidopsisSubNetwork_GeneList.txt Final_report.txt"
end


#--------------------------MAIN CODE-------------------------------------------------------------------------------

puts "Processing #{input_gene_list} file, this might take a while..."
Members.read_from_file(input_gene_list)


# Parameters for this assignment
INTACT_BASE_ADDRESS = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/'
SPECIES = 'species:arabidopsis' 
TAB25_FORMAT = 'tab25'
DEPTH = 2


Members.all_coexpresed_members.each do |gene|
  puts
  puts "Analyzing new gene #{gene.gene_id} with Uniprot ID #{gene.uniprot_id}"
  
  network = Networks.new  # create a seed network
  network.add_member(gene)
  network.recursive_search([gene], DEPTH) # search and assign all the interactors found to this net

  network.merge_with_common

  puts "This member is part of network #{network} with #{network.network_members.length} miembros"

end

Networks.reduce_networks


File.open(output_report_file, 'w') do |file|

  file.puts "This code was created by Miguel La Iglesia Mirones and Lucía Muñoz Gil"
  file.puts "Bioinformatics Programming Challenges Course at MSc Computational Biology (UPM)"
  file.puts "November 2023"
  file.puts "---------------------------------------------------------------------------------------"
  file.puts
  file.puts "Protein - protein interaction networks (interactome) final report for genes detailed" 
  file.puts "at #{input_gene_list} file"
  file.puts
  file.puts "---------------------------------------------------------------------------------------"
  file.puts
  file.puts "GLOBAL REPORT FOR ALL PROTEIN-PROTEIN INTERACTION NETWORKS"
  file.puts
  file.puts "Total number of nodes: #{Members.all_members.length}"
  file.puts "Number of nets: #{Networks.all_networks.length}"
  file.puts "Genes from list not included in any network: #{}"
  file.puts
  file.puts "---------------------------------------------------------------------------------------"
  file.puts
  file.puts "FEATURES OF EVERY INTERACTION NETWORK"
  Networks.all_networks.each_with_index do |network, idx|
    file.puts
    file.puts "Network ID #{idx + 1} (#{network}):"
    #file.puts "Network ID #{idx + 1} (#{network}) with members :"
    file.puts "  Number of members: #{network.network_members.length}"
    file.puts "  Genes from file part of this network:"
    
    network.network_members.each do |member|
      if member.gene_id
        file.puts "    Gene ID: #{member.gene_id}, Uniprot ID for coded protein: #{member.uniprot_id}" 
      end
    end

    file.puts
    file.puts "  KEGG annotations in network #{idx + 1}:"
    network.network_KEGGs.each do |key_, value|
      file.puts "    KEGG ID: #{key_}   Pathway: #{value}"
    end

    file.puts
    file.puts "  GO annotations in network #{idx + 1}:"
    network.network_GOs.each do |key_, value|
      file.puts "    GO ID: #{key_}   Term: #{value}"
    end

  end 
end


