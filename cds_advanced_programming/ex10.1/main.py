import example

print('First instance:')
instance_1 = example.Threading()
instance_1.run()
print('\n')

print('Second instance:')
instance_2 = example.Threading(3)
instance_2.run()

