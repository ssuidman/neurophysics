import numpy as np
import time
N = 10**9
t1 = time.time()
random_numbers = np.random.normal(0,1,N)
t2 = time.time()
print(round(t2-t1,4))
#running time is 25 sec