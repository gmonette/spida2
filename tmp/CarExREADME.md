# carEx


2019_01_22: copied spline and related functions and Rd files from spida2 --Georges

The main functions are:

- gsp: generates a spline function
- sc: generates a portion of a hypothesis (L) matrix to 
  estimate functions of spline parameters
- wald, walddf: estimate components and 
  perform simultaneous
  tests with L matrices even if rows are not 
  linearly independent
- Lfx: facilitates construction of L matrices
  for derivatives of continuous variables
- rowdiffs: pairwise differences of rows of
  L matrices to facilitate estimation of 
  differences between levels of categorical factors
- There are many functions to generate various 
  L matrices, most of which are probably obsolete
  
The most extensive examples are in the help files for 

- gsp
- wald, walddf
- Lfx
