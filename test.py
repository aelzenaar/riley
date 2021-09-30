import matlab.engine

matlab_eng = matlab.engine.start_matlab()
poly = [4,1]
print(matlab_eng.roots(poly))
