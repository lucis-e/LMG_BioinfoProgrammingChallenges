def recursive_search(depth, array_)

    if depth < 1
        return "end" 
    else
        puts "depth = #{depth}"
        array_.each do |element|
            puts element 
            modified_array = array_.map{|element| element + 1}
            puts "new array #{modified_array}"
        return recursive_search(depth-1, modified_array)
        end 
    end

end

list_of_nums=[4]


def recursive_search(depth, list_of_members)

    case depth
    when 0
        break

    else
        puts "#{depth}"
  
  end
