import sympy
from einsteinpy.symbolic import MetricTensor, ChristoffelSymbols, RiemannCurvatureTensor, RicciTensor, RicciScalar
sympy.init_printing()  # enables the best printing available in an environment
from sympy import simplify

#fill in metric
t,r,theta,phi = sympy.symbols('t r θ φ') #define symbols
nu = sympy.Function('ν')(t,r)
mu = sympy.Function('μ')(r)
metric_list = [[0 for i in range(4)] for i in range(4)] #create empty matrix
metric_list[0][0] = -sympy.exp(nu) #fill in the matrix (of the metric) with the symbols
metric_list[1][1] = sympy.exp(mu)
metric_list[2][2] = r**2
metric_list[3][3] = r**2*sympy.sin(theta)**2
# creating metric object
metric = MetricTensor(metric_list,[t,r,theta,phi]) #create metric from matrix (list-type object) and define what symbols are actually variables (so dt,dx,dy,dz)
metric.tensor() #metric.tensor() to show the tensor. metric.inv() to get inverse metric and metric.inv().tensor() to show inverse metric as tensor

#calculating other tensors starting with metric 
christoffel_symbols = ChristoffelSymbols.from_metric(metric)
riemann_tensor = RiemannCurvatureTensor.from_christoffels(christoffel_symbols) #=R[pho][sigma][mu][nu]
ricci_tensor = RicciTensor.from_riemann(riemann_tensor)
ricci_scalar = RicciScalar.from_riccitensor(ricci_tensor) #can also calculate it directly from metric, but if you want all objects then this is faster.

#To simplify expression you can do
#simplify(christoffel_symbols.tensor())



#Example runnning
#from einsteinpy.plotting import GeodesicPlotter
#from einsteinpy.examples import perihelion
#a = GeodesicPlotter()
#a.plot(perihelion())
#a.show()