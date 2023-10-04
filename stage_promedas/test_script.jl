# Create a large array of numbers
data = rand(10^8)

# Calculate the sum of squares using a loop
result = 0.0
t1 = time()
for x in data
    global result += x^2
end
t2 = time()

print(t2-t1)
