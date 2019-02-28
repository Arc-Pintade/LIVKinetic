rm nRoot.h
touch nRoot.h
cat nRoot.h >> ../include/kineticAnalyze.hpp
cat nRoot.h >> ../src/kineticAnalyze.cpp

#root -l -q TTstat.C

#vim ../stats/statsgTT_results.table
