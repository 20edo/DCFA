# DCFA
 
DCFA swept wing assignement

Teamwork
Team members: 

* Venti Edoardo         944421 
* Zemello Matteo        942003 
* Zucchelli Umberto     952952 

This simple program was developed during the 2019-2020 (during coronavirus pandemy) A.Y. as an assignement for the course of Dynamics and Control of Flexible Aircraft held by Pierangelo Masarati. 
The contributors are three students.
The scope of the assignement is to perform some aeroelastic analysis of a swept wing aircraft, our goal is to write a suite that can be easily extended to other engineering problems.

Brief introduction to the code

The code was originally developed to be as general as possible, though through the analysis some features have been built ad hoc. It is based on object-oriented programming.
The folder has been organized as follows: 
* the init script adds to path the folders containing the functions available for each object;
* Log folder contains the symbolic matlab scripts used to obtain matrices and vectors useful for the analysis. As a general rule they work in the element reference;
* Test folder contains some tests that may be used as examples to get acquainted with the code;
* Aircraft folder contains the model of the C-17 and the analysis performed;
* all the other folders have the name of an object and contain all the functions that work on that object and the object\textunderscore init function used to initialize the matlab structure associated with that object. The object\textunderscore init can be used as a reference to understand the meaning of the fields of the object structure.
\end{itemize}
