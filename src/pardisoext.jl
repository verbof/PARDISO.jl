# External functions to handle PARDISO functions

function pardisoSolveSystem(pardiso::ParDiSO, A::SparsePardisoCSR, b::Array)
    
    smbfctPARDISO(pardiso, A);
    factorPARDISO(pardiso, A);

    return solvePARDISO(pardiso, A, b);

end



