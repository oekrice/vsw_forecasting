#MPI (ish) CMA test. As the standard sweep approach is terrible, but this might be a decent compromise. Can have quite big parameter sets.

#Can work using

import cma
import numpy as np

popsize = 32
es = cma.CMAEvolutionStrategy(12 * [0], 0.5, {'verb_disp': 1, 'popsize': popsize})
while not es.stop():
   parameters = es.ask()
   #print(np.shape(parameters))
   solutions = [cma.ff.rosen(x) for x in parameters]
   #print(np.shape(solutions))
   print(np.min(solutions))
   es.tell(parameters, solutions)
