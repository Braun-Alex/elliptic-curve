package ec

import (
	"math/big"
	"strings"
)

// Parameters (XInSecp256k1G, YInSecp256k1G) for G and P in secp256k1

var XInSecp256k1G, _ = new(big.Int).SetString(
	strings.ToLower("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798"), 16)

var YInSecp256k1G, _ = new(big.Int).SetString(
	strings.ToLower("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8"), 16)

var P, _ = new(big.Int).SetString(
	strings.ToLower("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F"), 16)

type ElCPoint struct {
	X *big.Int
	Y *big.Int
}

// Returning base point G on the elliptic curve y^2 = x^3 + 7 (mod p) (secp256k1)

func BasePointGGet() (point ElCPoint) {
	return ElCPoint{XInSecp256k1G, YInSecp256k1G}
}

// Returning ElCPoint structure wrapped in coordinates

func ElCPointGen(x, y *big.Int) (point ElCPoint) {
	return ElCPoint{x, y}
}

// Checking that point is on curve y^2 = x^3 + 7 (mod p) (secp256k1)

func IsOnCurveCheck(a ElCPoint) (c bool) {
	XExpression := new(big.Int).Mul(a.X, a.X)
	XExpression.Mul(XExpression, a.X)
	YExpression := new(big.Int).Mul(a.Y, a.Y)
	result := new(big.Int).Sub(YExpression, XExpression)
	result.Sub(result, big.NewInt(7))
	result.Mod(result, P)
	return result.Cmp(big.NewInt(0)) == 0
}
