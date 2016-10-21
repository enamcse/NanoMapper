# NanoMapper
A Mapping tool using gapped minimizer and BWT FM-index

#Instructions:

- To compile the current version, give the following command from the current directory:
    
        g++ -o NanoMapper -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib ./Minimizer\ FM-index\ Merged\ Approach/Main_minimizer.cpp -lz -lsdsl -ldivsufsort -ldivsufsort64

- To run a demo, give the following command:

	    ./NanoMapper './DataSet/synthetic_reference.fasta' './DataSet/synthetic_reads.fasta' 14 24 './Results/Merged_enhanced_res_with_hash.txt' >> './Results/Summary.txt'

- To debug with short sequence, give the following command:

	    ./NanoMapper './DataSet/synthetic_reference_sample.fasta' './DataSet/synthetic_reads_sample.fasta' 14 24 './Results/Merged_enhanced_res_with_hash_sample.txt'  >> './Results/Summary.txt'

##Developer
Enamul Hassan and Md. Khairullah Gaurab under the direction of Ruhul Amin Shajib Sir

##Contact Information
Enamul Hassan : enamsustcse@gmail.com

Md. Khairullah Gaurab : mkgaurabsarkar@gmail.com

