# NanoMapper
A Mapping tool using gapped minimizer and BWT FM-index

#Prerequisites:

##API:
We have used [SDSL-lite](https://github.com/simongog/sdsl-lite) for fm-index. So, you need to install that first. To install that, you have to do the following:

- Make sure you have an updated [cmake](http://www.cmake.org/) build system.
- Make sure your `g++` version >= 4.9 and `clang` version >= 3.2
- Give the following commands in your terminal:
		
		git clone https://github.com/simongog/sdsl-lite.git
		cd sdsl-lite
		./install.sh

##OS:
- 64-bit unix based operating system (Linux or Mac OS X).

#Instructions:
- To get NanoMapper in your PC:

		git clone https://github.com/enamcse/NanoMapper.git && cd NanoMapper
		
- To compile the current version, give the following command from the current directory:
    
        g++ -o NanoMapper -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib ./Minimizer\ FM-index\ Merged\ Approach/Main_minimizer.cpp -lz -lsdsl -ldivsufsort -ldivsufsort64

- To run a demo, give the following command:

	    ./NanoMapper './DataSet/synthetic_reference.fasta' './DataSet/synthetic_reads.fasta' 14 24 2 './Results/Merged_enhanced_res.txt' >> './Results/Summary.txt'

##Developer
Enamul Hassan and Md. Khairullah Gaurab under the direction of Ruhul Amin Shajib Sir

##Contact Information
Enamul Hassan : enamsustcse@gmail.com
Md. Khairullah Gaurab : mkgaurabsarkar@gmail.com

