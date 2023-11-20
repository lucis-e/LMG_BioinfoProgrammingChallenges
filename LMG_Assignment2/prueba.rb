

prueba = [[0,1], [0,0]].select do |lista|


    lista.any? { |member| [9, 7, 3].include?(member)}  # check if there is any net with common members with the just created net
end

puts prueba.empty?

