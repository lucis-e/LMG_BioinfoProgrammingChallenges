require 'rest-client'
require 'json'

def togo_search(database,query,field="")
    togo_address = "http://togows.dbcls.jp/entry/#{database}/#{query}#{field}.json"
    puts togo_address
    togo_response = rest_api_request(togo_address)
    result = JSON.parse(togo_response.body)
    return result
end



def rest_api_request(address)
RestClient::Request.execute(
    method: :get,
    url: address,
    headers: {'Accept' => 'application/json'})  # use the RestClient::Request object's method "execute"
end

result = togo_search("uniprot", "AT2g13360","/accessions")

if result.is_a?(Array) && result.any?
    uniprot_id = result.first.first 
end 

puts uniprot_id