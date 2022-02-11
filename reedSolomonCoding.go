package main

import (
	"fmt"

	//"github.com/parksworks/gfGen.git"	// will be available..
	"myapp.com/gfGen" // local path will be replaced @ go.mod
)

type ReedSolomonCode struct {
	gf      gfGen.GaloisField
	genPoly gfGen.GfPoly

	n uint // code length
	k uint // message length
	p uint // number of parity symbols
}

// Make Generator polynomial of RS code over GF(2^m) with primitive polynomial
// [Inputs]
// n		Number of codeword symbols (code length)
// k		Number of message symbols
// m 		order of target galois field GF(2^m)
// priPoly	Primitive polynomial to generate GF(2^m)
// [Results] Make a generator polynomial of the RS code over GF(2^m) with the given primitive polynomial
func (rs *ReedSolomonCode) ConstructRSCode(n, k, m uint, priPoly *gfGen.GfPoly) error {
	// Constraints check
	if n <= k {
		return fmt.Errorf("[ReedSolomonCode] Cannot Construct RS Code! n should be larger than k!")
	} else if m <= 0 {
		return fmt.Errorf("[ReedSolomonCode] Cannot Construct RS Code! m should be a positive integer!")
	}
	// Construct Galois field
	gfErr := rs.gf.GfConstructGaloisField(m, priPoly)
	if gfErr != nil {
		return gfErr
	}

	// parameter setting and initialize
	rs.n = n
	rs.k = k
	rs.p = n - k // number of parity symbols
	if rs.genPoly == nil {
		rs.genPoly = make(gfGen.GfPoly)
	} else {
		for k2 := range rs.genPoly {
			delete(rs.genPoly, k2)
		}
	}

	// recursively multiply (x-a^i), for i = 1 ~ (n-k)
	rs.genPoly[0] = 0 // at first, 1*x^0
	for i := uint(1); i <= rs.p; i++ {
		tempPoly := rs.polyMulti(&rs.genPoly, &(gfGen.GfPoly{0: i, 1: 0})) // multiply by (a^i + x)
		// reset genPoly
		for k2 := range rs.genPoly {
			delete(rs.genPoly, k2)
		}
		// copy from tempPoly to genPoly
		for deg := range *tempPoly {
			rs.genPoly[deg] = (*tempPoly)[deg]
		}
	}

	return nil
}

// Do Encoding for a given data of length k, and return codeword of length n
// [Inputs]
// msg 		 		this slice should have k data in 0 ~ (2^m)-1 range
// [Results]		encoded RS codeword symbols of length n over GF(2^m)
func (rs *ReedSolomonCode) Encoding(msg *[]uint) (*[]uint, error) {
	if rs.genPoly == nil {
		return nil, fmt.Errorf("[ReedSolomonCode] Generator polynomial required! Please set parameters by using ConstructRSCode() method!")
	} else {
		msgPoly := rs.int2GfPoly(msg, rs.k)                // convert message to polynomial
		codewordPoly := rs.polyMulti(msgPoly, &rs.genPoly) // encoding using generator polynomial
		return rs.poly2Int(codewordPoly, rs.n), nil
	}
}

func (rs *ReedSolomonCode) Decoding(codeword *[]int) (*[]uint, error) {
	// check that the feasibility of decoding
	ersNum := uint(0)                       // number of erasures
	var ersPos []uint                       // position of erasures
	for deg := uint(0); deg < rs.n; deg++ { // deg should not be a negative int!
		if (*codeword)[deg] == -1 { // count the number of erasures
			ersPos = append(ersPos, uint(deg)) // append the position of erasures
			ersNum++
		} else if (*codeword)[deg] >= int(rs.gf.GfGetFieldSize()) {
			return nil, fmt.Errorf("[ReedSolomonCode] %dth codeword symbol %d is not a member of GF(%d)!", deg, (*codeword)[deg], rs.gf.GfGetFieldSize())
		}
	}
	if len(*codeword) > int(rs.n) { // length check
		return nil, fmt.Errorf("[ReedSolomonCode] Code length %d is not the same as %d! ", len(*codeword), rs.n)
	} else if ersNum > rs.p {
		return nil, fmt.Errorf("[ReedSolomonCode] Too many symbols (%d) are erased! This code can correct upto %d erasures!", ersNum, rs.p)
	}

	// get syndrome vector
	syndrome := make([]int, ersNum)
	for eqn := uint(0); eqn < ersNum; eqn++ { // it requires ersNum equations (is the same as the (power of root) - 1 of generator polynomial)
		syndrome[eqn] = -1
		for deg := uint(0); deg < rs.n; deg++ { // length of codeword is n
			if ((*codeword)[deg] != -1) && ((*codeword)[deg] != 0) { // -1: erased, 0: zero-coefficient (non-exponential form)4
				syndrome[eqn] = rs.addOverGf(syndrome[eqn], (*codeword)[deg]-1+int((eqn+1)*(deg))) // a^syndrome[eqn] + a^(codeword[deg]-1)*(a^(eqn+1)^deg) = a^symdrome[eqn] + a^(codeword[deg]-1+((eqn+1)*deg))
			}
		}
	}

	// get coefficient matrix
	coeffMatrix := make([][]int, ersNum)
	for eqn := uint(0); eqn < ersNum; eqn++ { // equation number (=power of primitive element that is a root of generator polynomial)
		coeffMatrix[eqn] = make([]int, ersNum)
		for ersInd := uint(0); ersInd < ersNum; ersInd++ { // erausre position index
			fmt.Printf("eqn = %d\tersInd = %d\t", eqn, ersInd)
			coeffMatrix[eqn][ersInd] = int((eqn+1)*ersPos[ersInd]) % (int(rs.gf.GfGetFieldSize()) - 1)
			fmt.Printf("%d\t", coeffMatrix[eqn][ersInd])
		}
		fmt.Println("")
	}

	fmt.Println("CoeffMat = ")
	for i := 0; i < int(ersNum); i++ {
		for j := 0; j < int(ersNum); j++ {
			fmt.Printf("%d\t", coeffMatrix[i][j])
		}
		fmt.Println("")
	}
	fmt.Print("Syndrome = ")
	for i := 0; i < int(ersNum); i++ {

		fmt.Printf("%d\t", syndrome[i])

	}
	fmt.Println("")
	// Gaussian elimination equation solving
	rs.gaussianElimination(&coeffMatrix, &syndrome)

	// construct decoding result
	decoded := make([]uint, rs.n)
	ersCnt := 0
	for deg := uint(0); deg < rs.n; deg++ {
		if (*codeword)[deg] != -1 {
			decoded[deg] = uint((*codeword)[deg])
		} else {
			decoded[deg] = uint(syndrome[ersCnt]) + 1
			ersCnt++
		}
	}
	return &decoded, nil
}

func (rs *ReedSolomonCode) addOverGf(a1, a2 int) int {
	if a1 == -1 {
		if a2 == -1 {
			return -1
		} else {
			return a2 % (int(rs.gf.GfGetFieldSize()) - 1)
		}
	} else if a2 == -1 {
		return a1 % (int(rs.gf.GfGetFieldSize()) - 1)
	} else {
		result := rs.gf.GfAdditionOfTwoExponents(uint(a1%(int(rs.gf.GfGetFieldSize())-1)), uint(a2%(int(rs.gf.GfGetFieldSize())-1)))
		fmt.Printf("ADDITION: a^%d + a^%d = a^%d\n", a1, a2, result)
		return result
	}
}

// multiply polynomials p1, p2
func (rs *ReedSolomonCode) polyMulti(p1, p2 *gfGen.GfPoly) *gfGen.GfPoly {
	var addedCoeff int
	fieldSize := rs.gf.GfGetFieldSize()
	returnPoly := make(gfGen.GfPoly)
	// p1 * p2
	for deg1, coef1 := range *p1 {
		for deg2, coef2 := range *p2 {
			if _, exist := returnPoly[deg1+deg2]; exist {
				addedCoeff = rs.gf.GfAdditionOfTwoExponents((coef1+coef2)%(fieldSize-1), returnPoly[deg1+deg2])
				if addedCoeff == -1 {
					delete(returnPoly, deg1+deg2)
				} else {
					returnPoly[deg1+deg2] = uint(addedCoeff)
				}
			} else {
				returnPoly[deg1+deg2] = (coef1 + coef2) % (fieldSize - 1)
			}
		}
	}
	return &returnPoly
}

// integer slice -> polynomial over Galois field
func (rs *ReedSolomonCode) int2GfPoly(data *[]uint, length uint) *gfGen.GfPoly {
	returnPoly := make(gfGen.GfPoly)
	for i := uint(0); i < length; i++ {
		if i != 0 {
			returnPoly[i] = i - 1
		}
	}
	return &returnPoly
}

// polynomial over Galois field -> integer slice
func (rs *ReedSolomonCode) poly2Int(poly *gfGen.GfPoly, length uint) *[]uint {
	returnInt := make([]uint, length)

	for deg, coeff := range *poly {
		returnInt[deg] = coeff + 1
	}

	return &returnInt
}

// Gaussian elimination equation solving
// This function assume that the coeffMatrix is a full-rank matrix! (based on the properties of RS codes)
// [Inputs]
// coeffMatrix *[][]uint 		coefficient matrix to be Gaussian eliminated
// syndrome *[]int				syndrome vector
// [Result]
// "syndrome" is directly modified, and it is the result
// [CAUTION]
// !!this function directly access and edit to memory of coeffMatrix and syndrome by pointer!!
func (rs *ReedSolomonCode) gaussianElimination(coeffMatrix *[][]int, syndrome *[]int) {
	syndromeNum := len(*syndrome)
	for row := 0; row < syndromeNum; row++ { // for each row
		invExp := (rs.gf.GfGetFieldSize() - 1) - uint((*coeffMatrix)[row][row]) // get inverse of syndrome value
		for col := row; col < syndromeNum; col++ {                              // make a^0 the leading element of the row
			(*coeffMatrix)[row][col] = ((*coeffMatrix)[row][col] + int(invExp)) % int(rs.gf.GfGetFieldSize()-1) // multiply inverse
		}
		if (*syndrome)[row] != -1 {
			(*syndrome)[row] = ((*syndrome)[row] + int(invExp)) % int(rs.gf.GfGetFieldSize()-1) // multiply the "inverse" corresponding syndrome vector
		}

		// eliminating the other rows
		for rowEliminator := 0; rowEliminator < syndromeNum; rowEliminator++ {
			if rowEliminator != row {
				multiplier := (*coeffMatrix)[rowEliminator][row] // get multiplier
				fmt.Printf("Selected multiplier = %d\n", multiplier)
				if multiplier != -1 {
					(*coeffMatrix)[rowEliminator][row] = -1
					for colEliminator := row + 1; colEliminator < syndromeNum; colEliminator++ {
						(*coeffMatrix)[rowEliminator][colEliminator] = rs.addOverGf((*coeffMatrix)[rowEliminator][colEliminator], (*coeffMatrix)[row][colEliminator]+int(multiplier))
					}
					if (*syndrome)[row] != -1 {
						(*syndrome)[rowEliminator] = rs.addOverGf((*syndrome)[rowEliminator], (*syndrome)[row]+int(multiplier))
					}
				}
			}
		}

		fmt.Println("CoeffMat = ")
		for i := 0; i < int(syndromeNum); i++ {
			for j := 0; j < int(syndromeNum); j++ {
				fmt.Printf("%d\t", (*coeffMatrix)[i][j])
			}
			fmt.Println("")
		}
		fmt.Print("Syndrome = ")
		for i := 0; i < int(syndromeNum); i++ {

			fmt.Printf("%d\t", (*syndrome)[i])

		}
		fmt.Println("")
	}
}

func main() {
	priPoly := &gfGen.GfPoly{0: 1, 1: 1, 3: 1} // 1x^0 + 1x^1 + 1x^3

	n := uint(7)
	k := uint(3)
	m := uint(3)

	var testRs ReedSolomonCode

	testRs.ConstructRSCode(n, k, m, priPoly)
	fmt.Printf("Generator Polynomial : ")
	for deg, coef := range testRs.genPoly {
		fmt.Printf("(a^%d * x^%d) + ", coef, deg)
	}

	msg := []uint{0, 1, 2}
	codeword, err := testRs.Encoding(&msg)

	fmt.Println("codeword")
	if err == nil {
		for i := 0; i < int(n); i++ {
			fmt.Printf("%d\t%d\n", i, (*codeword)[i])
		}
	}

	testRs.gf.GfPrintAddTable()

	// put erasures
	erased := make([]int, n)
	for i := 0; i < int(n); i++ {
		erased[i] = int((*codeword)[i])
	}
	erased[1] = -1
	erased[3] = -1
	erased[4] = -1
	erased[5] = -1

	fmt.Println("erased")
	for i := 0; i < int(n); i++ {
		fmt.Printf("%d\t%d\n", i, erased[i])
	}

	decoded, err2 := testRs.Decoding(&erased)

	if err2 != nil {
		fmt.Printf(err2.Error())
	} else {
		fmt.Println("index\tencocded\terased\tdecoded")
		for i := 0; i < int(n); i++ {
			fmt.Printf("%d\t%d\t%d\t%d\n", i, (*codeword)[i], erased[i], (*decoded)[i])
		}
	}

}
