package test

import (
	"Elliptic_Curve/pkg/ec"
	"crypto/rand"
	"math/big"
	"strings"
	"testing"
)

func SetRandom(bits int) *big.Int {
	randomNumber, err := rand.Int(rand.Reader, new(big.Int).Exp(big.NewInt(2),
		big.NewInt(int64(bits)), nil))
	if err != nil {
		panic("Could not be generated random " + string(rune(bits)) + "-bit number")
	}
	return randomNumber
}

func TestAssociativity(t *testing.T) {
	g := ec.BasePointGGet()
	k := SetRandom(256)
	d := SetRandom(256)
	h1 := ec.ScalarMult(*d, g)
	h2 := ec.ScalarMult(*k, h1)
	h3 := ec.ScalarMult(*k, g)
	h4 := ec.ScalarMult(*d, h3)
	if !ec.Eq(h2, h4) {
		t.Error("Group of elliptic curve points has no group associativity")
	}
}

func TestMainElCOperations(t *testing.T) {
	privateKey := new(big.Int)
	privateKey.SetString(strings.ToLower(
		"0D67243CBD68A96C95C849799F2B748CB641CE89A5D88C0652768272B11C2689"), 16)
	actualPublicKey := ec.ScalarMult(*privateKey, ec.BasePointGGet())
	expectedPublicKey := ec.StringToElCPoint(strings.ToLower(
		"0279BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798"))
	if !ec.Eq(actualPublicKey, expectedPublicKey) {
		t.Error("Elliptic curve functionality is not valid")
	}
}

func TestAdd(t *testing.T) {
	one := ec.ElCPoint{X: big.NewInt(8), Y: big.NewInt(3)}
	two := ec.ElCPoint{X: big.NewInt(3), Y: big.NewInt(6)}
	result := ec.AddElCPoints(one, two)
	expectedResult := ec.ElCPoint{X: big.NewInt(3), Y: big.NewInt(5)}
	if !ec.Eq(result, expectedResult) {
		t.Error(result, "\n", expectedResult)
	}
}

func TestMul(t *testing.T) {
	point := ec.ElCPoint{X: big.NewInt(1), Y: big.NewInt(8)}
	result := ec.DoubleElCPoints(point)
	expectedResult := ec.ElCPoint{X: big.NewInt(7), Y: big.NewInt(7)}
	if !ec.Eq(result, expectedResult) {
		t.Error(result, "\n", expectedResult)
	}
}

func TestScalarMult(t *testing.T) {
	point := ec.ElCPoint{X: big.NewInt(9), Y: big.NewInt(4)}
	result := ec.ScalarMult(*big.NewInt(5), point)
	expectedResult := ec.ElCPoint{X: big.NewInt(9), Y: nil}
	if !ec.Eq(result, expectedResult) {
		t.Error(result, "\n", expectedResult)
	}
}
