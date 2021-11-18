import sympy
from einsteinpy.symbolic import MetricTensor, ChristoffelSymbols, RiemannCurvatureTensor, RicciTensor, RicciScalar
sympy.init_printing()  # enables the best printing available in an environment
from sympy import simplify

#fill in metric
t,x,y,z,w = sympy.symbols('t x y z w') #define symbols
metric_list = [[0 for i in range(4)] for i in range(4)] #create empty matrix
metric_list[0][0] = -(1-w**2*(x**2+y**2)) #fill in the matrix (of the metric) with the symbols
metric_list[0][1] = -2*w*y
metric_list[1][0] = -2*w*y
metric_list[0][2] = 2*w*x
metric_list[2][0] = 2*w*x
metric_list[1][1] = 1
metric_list[2][2] = 1
metric_list[3][3] = 1
# creating metric object
metric = MetricTensor(metric_list,[t,x,y,z]) #create metric from matrix (list-type object) and define what symbols are actually variables (so dt,dx,dy,dz)
metric.tensor() #metric.tensor() to show the tensor. metric.inv() to get inverse metric and metric.inv().tensor() to show inverse metric as tensor

#calculating other tensors starting with metric 
christoffel_symbols = ChristoffelSymbols.from_metric(metric)
riemann_tensor = RiemannCurvatureTensor.from_christoffels(christoffel_symbols)
ricci_tensor = RicciTensor.from_riemann(riemann_tensor)
ricci_scalar = RicciScalar.from_riccitensor(ricci_tensor) #can also calculate it directly from metric, but if you want all objects then this is faster.

#To simplify expression you can do
#simplify(christoffel_symbols[2])
