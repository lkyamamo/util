#FLAGS=-D NNQMD
FLAGS=-g -std=c++17 -fsanitize=address 
#FLAGS=-g -std=c++17 
#FLAGS=-g -std=c++17 

corr: corr.cpp
	g++ ${FLAGS} -D NNQMD -I${BOOST_ROOT}/include corr.cpp -o nnqmd_corr.exe &
	g++ ${FLAGS} -I${BOOST_ROOT}/include corr.cpp -o qmd_corr.exe

msd: msd.cpp
	g++ -D NNQMD ${FLAGS} msd.cpp -o nnqmd_msd.exe &
	g++ ${FLAGS} msd.cpp -o qmd_msd.exe

h:
	g++ ${FLAGS} -I${BOOST_ROOT}/include histogram.cpp -o hist.exe

nbr:
	g++ ${FLAGS} nbr.cpp -o nbr.exe

