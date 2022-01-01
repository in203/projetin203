CC=mpicc
CXX=mpic++

CXXFLAGS+=$(mpic++ --showme:compile)
LIBS+=$(mpic++ --showme:link)
MPIOPT= -oversubscribe

