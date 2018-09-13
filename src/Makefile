
CXX = g++
CXXFLAGS += -Wall -std=c++11 -I$(mkEigenInc)

EXE = main

all: $(EXE)

.PHONY: all clean distclean

OBJS = Helpers.o Distance.o Coordinates.o Distance_Tplane.o EmpiricalVariogram.o FittedVariogram.o VariogramFactoryProxy.o DesignMatrix.o DistanceManifold.o MapFunctions.o Model.o

$(EXE) : $(OBJS) main.o
	$(CXX) $(CXXFLAGS) $^ -o $@



Distance.o: Distance.hpp
Coordinates.o: Coordinates.hpp Helpers.hpp Distance.hpp
Distance_Tplane.o: Distance_Tplane.hpp Helpers.hpp
Helpers.o: Helpers.hpp
VariogramFactoryProxy.o: HelpersFactory.hpp FittedVariogram.hpp
FittedVariogram.o: EmpiricalVariogram.hpp FittedVariogram.hpp
EmpiricalVariogram.o: EmpiricalVariogram.hpp Distance_Tplane.hpp Coordinates.hpp Distance.hpp Helpers.hpp
main2.o: Helpers.hpp FittedVariogram.hpp EmpiricalVariogram.hpp HelpersFactory.hpp
DistanceManifold.o: DistanceManifold.hpp Helpers.hpp
MapFunctions.o: MapFunctions.hpp Distance_Manifold.hpp
DesignMatrix.o: DesignMatrix.hpp Coordinates.hpp
Model.o: Model.hpp Helpers.hpp

#Point.o Coordinates.o Distance.o: Point.hpp
#Coordinates.o Variogram.o: Coordinates.hpp
#Distance.o Variogram.o: Distance.hpp
#Distance_Tplane.o Variogram.o: Distance_Tplane.hpp
#Variogram.o: Variogram.hpp
#Helpers.o Variogram.o Point.o Distance_Tplane.o Coordinates.o: Helpers.hpp
#Helpers.o: Proxy.hpp
#Helpers.o: Factory.hpp

clean:
	$(RM) *.o
	$(RM) $(EXE)
distclean: clean
	$(RM) *.~
