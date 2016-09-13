# External functions to handle PARDISO functions

function pardisoSolveSystem(pardiso::ParDiSo, A::SparsePardisoCSR, b::Array)
    
    smbfctPARDISO(pardiso, A);
    factorPARDISO(pardiso, A);

    return solvePARDISO(pardiso, A, b);

end


function pardisoFactor(mtype, A::SparsePardisoCSR)

    pardiso = ParDiSo(mytpe, 0);

    smbfctPARDISO(pardiso, A);
    factorPARDISO(pardiso, A);

    return pardiso

end


function pardisoSelInv(pardiso::ParDiSo, A::SparsePardisoCSR)
    
    smbfctPARDISO(pardiso, A);
    factorPARDISO(pardiso, A);

    return invertPARDISO(pardiso, A);

end


pardisoSelInv(pardiso::ParDiSo, A::SparseMatrixCSC) = pardisoSelInv(pardiso, SparsePardisoCSR(A))

