testthat::test_that(
    desc = "getMolecule handles z parameter (charge) correctly", 
    code = {
        
        suppressMessages({
            # exactmass is corrected for electron mass if z!=0 specified
            testthat::expect_equal(getMolecule("H", z=0)$exactmass, 1.00782503)
            testthat::expect_equal(getMolecule("H", z=1)$exactmass, 1.00727645)
            
            # also isotope masses need to be adjusted
            testthat::expect_equal(getMolecule("C6H12O6", z=2)$isotopes[[1]][1,1], 90.031146)
        })
        
    }
)

testthat::test_that(
    desc = "getMonoisotopic returns monoisotopic masses correctly", 
    code = {
        
        # example where exactmass and monoisotopic mass differ
        fml <- "C89H176O16P2"
        out <- getMolecule(fml)
        
        # the exact mass is the second isotope
        testthat::expect_equal(getMass(out), getIsotope(out, 2)[[1]][1,])
        
        # the monoisotopic mass is the first isotope
        testthat::expect_equal(unname(getMonoisotopic(out)), getIsotope(out, 1)[[1]][1,])
        
        # getMonoisotopic works on a fml directly
        testthat::expect_equal(getMonoisotopic(out), getMonoisotopic(fml))
        
        # getMonoisotopic works also on decomposeIsotopes output
        out <- decomposeIsotopes(c(147.0529, 148.0563), c(100.0, 5.56))
        testthat::expect_equal(unname(getMonoisotopic(out)), sapply(getIsotope(out, 1), function(x) { x[1,] }))
    }
)

testthat::test_that(
    desc = "getMolecule checks sum formulas to prevent crashes", 
    code = {
        
        # this was an error reported in issue #3
        # strings which can not be processed by the C++ function return an error
        testthat::expect_error(getMolecule("3"))
        
        # elements contained in the formula need to exist in the provided elements list
        testthat::expect_error(getMolecule(formula = "CH4", elements = initializeElements("H")))

    }
)