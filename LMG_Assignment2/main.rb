#----------------- LOAD EXTERNAL LIBRARIES AND MODULES-------------------------------------------------------------------
require 'rest-client'
require 'json'
require './Members'
require './Networks'
require './AnnotatedMembers'

#------------------------ PUBLIC INSTANCE METHODS -----------------------------------------------------------------------

# Programatic REST API get request given an http address
#
# @param address [String] http/url address to get request
# @return [void]
def rest_api_request(address)
  RestClient::Request.execute(
    method: :get,
    url: address,
    headers: {'Accept' => 'application/json'})
end

# Process tab25 column syntax structure to extract value: <XREF>:<VALUE>(<DESCRIPTION>). Extracts VALUE without quotes if any.
#
# @param element [String] selected tab-delimited column from REST API response tab25 file
# @return [String] VALUE without quotes if any
def extract_xref(element)
  match_data = element.match(/:(\S+)(?:\(|$)/)
  val = match_data[1].gsub(/\A"/, '').gsub(/"\z/, '') if match_data
  return(val)
end

# Get request for multiple databases of TOGO REST API recieved in json format and parsing.
#
# @param database [String] Name of the database
# @param query [String] Query (e.g. identifier) to get response from REST API
# @param field [String] Specific databases' directory to get response from
# @return [String] Parsed JSON object as a hash with the TOGO REST API response 
def togo_search(database,query,field="")
  togo_address = "http://togows.dbcls.jp/entry/#{database}/#{query}#{field}.json"
  togo_response = rest_api_request(togo_address)
  result = JSON.parse(togo_response.body)
  return result
end

#------------------------------ INPUT FILES NAMES AS ARGUMENTS FROM  COMMAND LINE ---------------------------------------

# Check for one common error: not specifying the correct number of files needed for the program to run
if ARGV.length != 2
  abort "Incorrect number of files passed. Two files names must be specified: input list of genes and name for final_report"
end

# Check for second common error: incorrect usage, files in incorrect order or wrong name passed.
if ARGV[0] == "ArabidopsisSubNetwork_GeneList.txt" && ARGV[1] == "Final_report.txt"
  input_gene_list = ARGV[0]
  output_report_file = ARGV[1]
else 
  abort "Incorrect order of files passed. Usage: ruby main.rb ArabidopsisSubNetwork_GeneList.txt Final_report.txt"
end

#-------------------------------------- SET PARAMETERS FOR THIS ASSIGNMENT ----------------------------------------------
# Intact DataBase base address
INTACT_BASE_ADDRESS = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/'

# Species to which all proteins found in the protein-protein interaction response should belong
SPECIES = 'species:arabidopsis' 

# Format for Intact database protein-protein interaction response
TAB25_FORMAT = 'tab25'

# Depth of recursive interactor search (how many protein-protein interaction searches from initial "coexpressed" gene)
DEPTH = 2

#------------------------------------------------- MAIN CODE ------------------------------------------------------------
puts "Processing #{input_gene_list} file, this might take a while..."

# Process file with "co-expressed" gene's locus names specified from command line. Create Member instances for all genes in file
LuMike_objects::Members.read_from_file(input_gene_list)

# Iterate through the array with just-created Member instances for the "coexpressed" genes
# Search for its direct interactors and the direct interactors of them
# All found interactors + the initial gene are part of the same network
# Finally, networks that are conected to other are actually the same net, so should be merged
LuMike_objects::Members.all_coexpresed_members.each do |gene|
  puts
  puts "Analyzing new gene #{gene.gene_id} with Uniprot ID #{gene.uniprot_id}"
  network = LuMike_objects::Networks.new  # Create new Networks class instance
  network.add_member(gene)  # Add coexpressed gene Member instance as member of the network
  network.recursive_search([gene], DEPTH)   # Protein-protein interactor recursive search for 2 levels of depth, found interactors would be added as part of the network
  network.merge_with_common   # Check for common members with other already-created networks and merge if possitive 
  puts "This member is part of network #{network} with #{network.network_members.length} member(s)"
end

# Delete networks with only one member (the "coexpressed" gene) since one node it is not a net
LuMike_objects::Networks.reduce_networks

# Get KEGG and GO annotations for all members of network and asign its UNION to network
LuMike_objects::Networks.all_networks.each do |network|
  network.annotate_network
end

# Create final report with findings of out analysis
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
  file.puts "Total number of nodes: #{LuMike_objects::Members.number_of_members}"
  file.puts "Number of nets: #{LuMike_objects::Networks.get_number_of_nets}"
  file.puts "Genes from list not included in any network: #{LuMike_objects::Networks.nodes_without_net}"
  file.puts
  file.puts "---------------------------------------------------------------------------------------"
  file.puts
  file.puts "FEATURES OF EVERY INTERACTION NETWORK"
  LuMike_objects::Networks.all_networks.each_with_index do |network, idx|
    file.puts
    file.puts "Network ID #{idx + 1} (#{network}):"
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