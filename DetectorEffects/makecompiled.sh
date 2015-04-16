g++ `root-config --cflags` src/DetectorEffects.C -c
g++ `root-config --cflags` MyfarmoutAnalyzer.cc DetectorEffects.o -o main.exe `root-config --libs`
#g++ `root-config --cflags` Main.C DetectorEffects.o -o main.exe `root-config --libs`
