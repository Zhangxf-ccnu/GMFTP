README file for Matlab code supporting the paper "Detecting overlapping protein complexes based on a generative model of functional and topological properties".


Contents of this archive
------------------------
This archive contains several Matlab scripts used to detect protein complexes using the algorithm GMFTP described in the above paper. 

(1) GMFTP.m: Matlab script for the core of GMFTP which detects protein complexes from a given PPI network with adjacent matrix A and a functional profile with association matrix F.
 
(2) GMFTP_main.m: Matlab script for the main function of GMFTP which read data from text files first and then detect protein complexes using GMFTP.m. It also writes the detected complexes into a file "output_file_name" provided by the user or the default file "complex_result.txt"

(3) GMFTP_demo.m: A simple Matlab script to test GMFTP. When a data set is chose, it can be run in a straightforward manner within a Matlab window.

This archive also contains a folder named as "data" which includes the six PPI networks and the corresponding total GO annotations used in this study.

Please note that due to the local minimum with different random initializations, the result may change a little. Therefore, to guard against the possibility of getting stuck in a local 
minimum, we suggest repeating the entire calculations multiple times and choosing the result that gives the lowest value of objective function.

Please do not hesitate to contact Prof. Dao-Qing Dai at stsddq@mail.sysu.edu.cn (or Dr. Xiao-Fei Zhang at zhangxf9@mail2.sysu.edu.cn) to seek any clarifications regarding any contents or operation of the archive.
