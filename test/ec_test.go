package test

import (
	"elliptic-curve/pkg/ec"
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

func TestGroupOperation(t *testing.T) {
	g := ec.BasePointGGet()
	randomFirstScalar := SetRandom(256)
	randomSecondScalar := SetRandom(256)
	randomFirstPoint := ec.ScalarMult(*randomFirstScalar, g)
	randomSecondPoint := ec.ScalarMult(*randomSecondScalar, g)
	resultPoint := ec.AddElCPoints(randomFirstPoint, randomSecondPoint)
	if !ec.IsOnCurveCheck(resultPoint) {
		t.Error("Group operation result does not belong to the elliptic curve")
	}
}

func TestGroupAssociativity(t *testing.T) {
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

func TestGroupNeutralElement(t *testing.T) {
	g := ec.BasePointGGet()
	randomScalar := SetRandom(256)
	randomPoint := ec.ScalarMult(*randomScalar, g)
	xRandom := SetRandom(256)
	infinitePoint := ec.ElCPointGen(xRandom, nil)
	leftPoint := ec.AddElCPoints(randomPoint, infinitePoint)
	rightPoint := ec.AddElCPoints(infinitePoint, randomPoint)
	if !ec.Eq(leftPoint, rightPoint) || !ec.Eq(randomPoint, leftPoint) ||
		!ec.Eq(randomPoint, rightPoint) {
		t.Error("Group of elliptic curve points has no neutral element that a*0 = 0*a = a")
	}
}

func TestGroupInverseElement(t *testing.T) {
	g := ec.BasePointGGet()
	randomScalar := SetRandom(256)
	randomPoint := ec.ScalarMult(*randomScalar, g)
	inverseRandomPoint := ec.ElCPointGen(new(big.Int).Set(randomPoint.X),
		new(big.Int).Mul(randomPoint.Y, big.NewInt(-1)))
	resultPoint := ec.AddElCPoints(randomPoint, inverseRandomPoint)
	if resultPoint.Y != nil {
		t.Error("Group of elliptic curve points has no inverse element that a*a^(-1) = 0")
	}
}

func TestGroupCommutativity(t *testing.T) {
	g := ec.BasePointGGet()
	randomFirstScalar := SetRandom(256)
	randomSecondScalar := SetRandom(256)
	randomFirstPoint := ec.ScalarMult(*randomFirstScalar, g)
	randomSecondPoint := ec.ScalarMult(*randomSecondScalar, g)
	if !ec.Eq(ec.AddElCPoints(randomFirstPoint, randomSecondPoint),
		ec.AddElCPoints(randomSecondPoint, randomFirstPoint)) {
		t.Error("Group of elliptic curve points has no group commutativity")
	}
}

func TestMainElCOperations(t *testing.T) {
	privateKey := new(big.Int)
	privateKey.SetString(strings.ToLower(
		"0D67243CBD68A96C95C849799F2B748CB641CE89A5D88C0652768272B11C2689"), ec.HexEncoding)
	actualPublicKey := ec.ScalarMult(*privateKey, ec.BasePointGGet())
	expectedPublicKey := ec.StringToElCPoint(strings.ToLower(
		"0324BD45AE802ACCF2FCF7D9E2639356F0CE5C21704C26943D130D2708B26B80E6"))
	if !ec.Eq(actualPublicKey, expectedPublicKey) {
		t.Error("Elliptic curve operations have not been properly implemented")
	}
}

func TestIsElCPointOnTheCurve(t *testing.T) {
	x, _ := new(big.Int).SetString(strings.ToLower(
		"24BD45AE802ACCF2FCF7D9E2639356F0CE5C21704C26943D130D2708B26B80E6"), ec.HexEncoding)
	y, _ := new(big.Int).SetString(strings.ToLower(
		"ABA2A81A7616C5FF79CCA277214BDA6F1A04CF24145D9713EE7D0CD319C16A9B"), ec.HexEncoding)
	point := ec.ElCPointGen(x, y)
	if !ec.IsOnCurveCheck(point) {
		t.Error("Elliptic curve point is not on the elliptic curve")
	}
}

func TestEncodingOfElCPoint(t *testing.T) {
	x, _ := new(big.Int).SetString(strings.ToLower(
		"24BD45AE802ACCF2FCF7D9E2639356F0CE5C21704C26943D130D2708B26B80E6"), ec.HexEncoding)
	y, _ := new(big.Int).SetString(strings.ToLower(
		"ABA2A81A7616C5FF79CCA277214BDA6F1A04CF24145D9713EE7D0CD319C16A9B"), ec.HexEncoding)
	point := ec.ElCPointGen(x, y)
	encodedPoint := ec.ElCPointToString(point)
	if encodedPoint != strings.ToLower(
		"0324BD45AE802ACCF2FCF7D9E2639356F0CE5C21704C26943D130D2708B26B80E6") {
		t.Error("Encoded elliptic curve point has not been properly encoded")
	}
}

func TestDecodingOfElCPoint(t *testing.T) {
	decodedPoint := ec.StringToElCPoint(strings.ToLower(
		"0324BD45AE802ACCF2FCF7D9E2639356F0CE5C21704C26943D130D2708B26B80E6"))
	x, _ := new(big.Int).SetString(strings.ToLower(
		"24BD45AE802ACCF2FCF7D9E2639356F0CE5C21704C26943D130D2708B26B80E6"), ec.HexEncoding)
	y, _ := new(big.Int).SetString(strings.ToLower(
		"ABA2A81A7616C5FF79CCA277214BDA6F1A04CF24145D9713EE7D0CD319C16A9B"), ec.HexEncoding)
	actualPoint := ec.ElCPointGen(x, y)
	if !ec.Eq(decodedPoint, actualPoint) {
		t.Error("Encoded elliptic curve point has not been properly decoded")
	}
}

func TestConvertingOfElCPoint(t *testing.T) {
	x, _ := new(big.Int).SetString(strings.ToLower(
		"24BD45AE802ACCF2FCF7D9E2639356F0CE5C21704C26943D130D2708B26B80E6"), ec.HexEncoding)
	y, _ := new(big.Int).SetString(strings.ToLower(
		"ABA2A81A7616C5FF79CCA277214BDA6F1A04CF24145D9713EE7D0CD319C16A9B"), ec.HexEncoding)
	actualPoint := ec.ElCPointGen(x, y)
	encodedPoint := ec.ElCPointToString(actualPoint)
	decodedPoint := ec.StringToElCPoint(encodedPoint)
	if !ec.Eq(decodedPoint, actualPoint) {
		t.Error("Converting of elliptic curve point has not been properly implemented")
	}
}
