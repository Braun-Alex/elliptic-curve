package ec

import (
	"fmt"
	"math/big"
	"strings"
)

var OneInASCII byte = 49

var HexEncoding int = 16

// Elliptic curve y^2 = x^3 + 7 (mod p)

var A = big.NewInt(0)

var B = big.NewInt(7)

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

func BasePointGGet() ElCPoint {
	return ElCPoint{XInSecp256k1G, YInSecp256k1G}
}

// Returning ElCPoint structure wrapped in coordinates

func ElCPointGen(x, y *big.Int) ElCPoint {
	return ElCPoint{x, y}
}

// Checking that point is on curve y^2 = x^3 + 7 (mod p) (secp256k1)

func IsOnCurveCheck(a ElCPoint) bool {
	XExpression := new(big.Int).Mul(a.X, a.X)            // x^2 = x*x
	XExpression.Mul(XExpression, a.X)                    // x^3 = x^2*x
	YExpression := new(big.Int).Mul(a.Y, a.Y)            // y^2 = y*y
	result := new(big.Int).Sub(YExpression, XExpression) // result := y^2 - x^3
	result.Sub(result, new(big.Int).Mul(A, a.X))         // result -= a*x
	result.Sub(result, B)                                // result -= b
	result.Mod(result, P)                                // result %= p
	return result.Cmp(big.NewInt(0)) == 0
}

// Adding two different elliptic curve points

func AddElCPoints(a, b ElCPoint) ElCPoint {
	if a.X == nil || b.X == nil { // Elliptic curve points must be valid
		return ElCPoint{nil, nil}
	} else if a.Y == nil && b.Y == nil { // O(x, inf) + O(x, inf) = O(x, inf)
		averageX := new(big.Int).Add(a.X, b.X)
		averageX.Div(averageX, big.NewInt(2))
		return ElCPoint{averageX, nil}
	} else if a.Y == nil { // O(x, inf) + P(x, y) = P(x, y)
		return b
	} else if b.Y == nil { // P(x, y) + O(x, inf) = P(x, y)
		return a
	} else if a.X.Cmp(b.X) == 0 && a.Y.Cmp(b.Y) == 0 { // P(x, y) + P(x, y) = 2P(x, y)
		return DoubleElCPoints(a)
	} else if a.X.Cmp(b.X) == 0 && new(big.Int).Add(a.Y, b.Y).Cmp(big.NewInt(0)) == 0 {
		return ElCPoint{new(big.Int).Set(a.X), nil} // P(x, y) + P(x, -y) = O(x, inf)
	} else {
		difference := new(big.Int).Sub(b.Y, a.Y)               // d = y2 - y1
		difference.Div(difference, new(big.Int).Sub(b.X, a.X)) // d /= (x2 - x1)
		lambda := difference.Mod(difference, P)                // d %= p; λ = d
		result := ElCPoint{new(big.Int), new(big.Int)}         // r := ElCPoint{x3, y3}
		result.X = new(big.Int).Mul(lambda, lambda)            // λ^2 = λ*λ; x3 = λ
		result.X.Sub(result.X, a.X)                            // x3 -= x1
		result.X.Sub(result.X, b.X)                            // x3 -= x2
		result.X.Mod(result.X, P)                              // x3 %= p
		result.Y.Mul(lambda, new(big.Int).Sub(a.X, result.X))  // y3 = λ*(x1 - x3)
		result.Y.Sub(result.Y, a.Y)                            // y3 -= y1
		result.Y.Mod(result.Y, P)                              // y3 %= p
		return result
	}
}

// Double multiplying of elliptic curve point

func DoubleElCPoints(a ElCPoint) ElCPoint {
	if a.X == nil { // Elliptic curve point must be valid
		return ElCPoint{nil, nil}
	} else if a.Y == nil { // 2*O(x, inf) = O(x, inf)
		return ElCPoint{new(big.Int).Add(a.X, a.X), nil}
	} else {
		lambda := new(big.Int).Mul(a.X, a.X)                  // λ = x1^2
		lambda.Mul(lambda, big.NewInt(3))                     // λ *= 3
		lambda.Add(lambda, A)                                 // λ += a
		lambda.Div(lambda, new(big.Int).Add(a.Y, a.Y))        // λ /= (x2 - x1)
		lambda.Mod(lambda, P)                                 // λ %= p
		result := ElCPoint{new(big.Int), new(big.Int)}        // r := ElCPoint{x3, y3}
		result.X = new(big.Int).Mul(lambda, lambda)           // λ^2 = λ*λ; x3 = λ
		result.X.Sub(result.X, new(big.Int).Add(a.X, a.X))    // x3 -= 2*x1
		result.X.Mod(result.X, P)                             // x3 %= p
		result.Y.Mul(lambda, new(big.Int).Sub(a.X, result.X)) // y3 = λ*(x1 - x3)
		result.Y.Sub(result.Y, a.Y)                           // y3 -= y1
		result.Y.Mod(result.Y, P)                             // y3 %= p
		return result
	}
}

// Scalar multiplying of elliptic curve point

func ScalarMult(k big.Int, a ElCPoint) ElCPoint {
	var resultPoint ElCPoint
	if k.Cmp(big.NewInt(0)) == 0 { // 0*P(x, y) = O(x, inf)
		return ElCPoint{a.X, nil}
	} else if k.Cmp(big.NewInt(0)) == -1 { // -kP(x, y) = kP(x, -y)
		return ScalarMult(*new(big.Int).Mul(&k, big.NewInt(-1)),
			ElCPoint{a.X, new(big.Int).Mul(a.Y, big.NewInt(-1))})
	} else if k.Cmp(big.NewInt(1)) == 0 {
		return a
	} else {
		bits := k.Text(2)                           // Binary representation of decimal scalar
		doublePoints := make([]ElCPoint, len(bits)) // Allocation of slice with length of number of bits
		for _, point := range doublePoints {        // P(x, y) + O(x, inf) = P(x, y)
			point.X = big.NewInt(-1)
			point.Y = nil
		} // All double points in slice before computation are O(x, inf)
		bufferPoint := a
		bufferIndex := 0
		for i := len(bits) - 1; i >= 0; i-- {
			if bits[i] == OneInASCII { // Symbol bits[i] is "1"
				doublePoints[i] = bufferPoint
				for l := len(bits) - 1 - i; l > bufferIndex; l-- {
					doublePoints[i] = DoubleElCPoints(doublePoints[i])
				}
				bufferPoint = doublePoints[i]
				bufferIndex = len(bits) - 1 - i
			}
		}
		for i := range doublePoints { // Initializing result
			if doublePoints[i].Y != nil {
				resultPoint = doublePoints[i]
				bufferIndex = i
				break
			}
		} // If bit is 0, the correspondent double point is O(x, inf); R(x, y) + O(x, inf) = R(x, y)
		for i := bufferIndex + 1; i < len(doublePoints); i++ {
			resultPoint = AddElCPoints(resultPoint, doublePoints[i])
		} // Result will change while adding when the double point is not equal to O(x, inf)
	}
	return resultPoint
}

func ElCPointToString(point ElCPoint) string {
	var leastSignificantBit string // Indicates that number is even or odd
	if point.X == nil {
		return "Not valid point"
	}
	if point.Y == nil {
		return "Infinite point" // O(x, inf)
	}
	if point.Y.Bit(point.Y.BitLen()-1) == 0 {
		leastSignificantBit = "02" // 0
	} else {
		leastSignificantBit = "03" // 1
	}
	return leastSignificantBit + point.X.Text(HexEncoding)
}

func StringToElCPoint(s string) ElCPoint {
	leastSignificantBit := s[:2]
	if leastSignificantBit != "02" && s[:2] != "03" {
		return ElCPoint{nil, nil} // Not valid point
	}
	l := new(big.Int)
	if leastSignificantBit == "02" {
		l.SetInt64(0)
	} else {
		l.SetInt64(1)
	}
	x, y := new(big.Int), new(big.Int)
	x.SetString(s[2:], HexEncoding)
	y.Exp(x, big.NewInt(3), nil)     // x^3
	y.Add(y, new(big.Int).Mul(x, A)) // x^3 + a*x
	y.Add(y, B)                      // x^3 + a*x + b
	hasRoot := y.ModSqrt(y, P)       // Square root of x^3 + a*x + b modulo p
	if hasRoot == nil {
		return ElCPoint{nil, nil} // Such point is not on elliptic curve
	}
	if l.Cmp(new(big.Int).Mod(y, big.NewInt(2))) == 0 {
		return ElCPoint{x, y}
	}
	return ElCPoint{x, y.Sub(P, y)}
}

func PrintElCPoint(point ElCPoint) {
	fmt.Printf("Elliptic curve point has such coordinates: \n%s\n%s\n",
		point.X.Text(HexEncoding), point.Y.Text(HexEncoding))
}
