min = gets.to_f

puts "min = #{min}"

results = [min]

limit = min * 1.02

puts "limit = #{limit}"

n = 1.0

while x = gets
  results.push(x.to_f) if x.to_f <= limit
  n += 1.0
end

puts "These are the results below the limit:"
puts results

p = results.size / n

x = (1 - (1 - p) ** n) * 100

puts "n = #{n}"
puts "p = #{p}"
puts "x = #{x}"
