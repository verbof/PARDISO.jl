#=
n = 10;
m = 8;
rowptr = int32([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12]);
colptr = int32([1, 8, 7, 3, 8, 5, 5, 6, 1, 3, 7]);
val    = [float64(i) for i in 1:11];


A = SparseMatrixCSC(m, n, colptr, rowptr, val);


rowptr64 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12];
colptr64 = [1, 8, 7, 3, 8, 5, 5, 6, 1, 3, 7];


B = SparseMatrixCSC(m, n, colptr64, rowptr64, val);
=#

#=

     [  2.9  0   0.5 -0.2 0   ]
     [  0    1.0 0    0   0   ]
X =  [  0    0   2.0  0   0   ]
     [  0    0   0    1.2 0   ]
     [  0.9 -0.3 0.4  0   3.0 ]

CSR FORMAT:
    nzvalues: 2.9 0.5 -0.2 1.0 2.0 1.2 0.9 -0.3 0.4 3.0
    colind:   1   3    4   2   3   4   1    2   3   5
    
    rowptr:   1 4 5 6 7 11
***


     [  2.9  0   0    0   0.9 ]
     [  0    1.0 0    0  -0.3 ]
X' = [  0.5  0   2.0  0   0.4 ]
     [ -0.2  0   0    1.2 0   ]
     [  0    0   0    0   3.0 ]

CSC FORMAT:
    nzvalues: 2.9 0.5 -0.2 1.0 2.0 1.2 0.9 -0.3 0.4 3.0
    rowind:   1   3    4   2   3   4   1    2   3   5
    
    colptr:   1 4 5 6 7 11
***

=#

using PARDISO;

pardiso = ParDiSo(-2,1);

initPARDISO(pardiso);
