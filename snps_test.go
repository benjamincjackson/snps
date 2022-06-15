package main

import (
	"bytes"
	"fmt"
	"strings"
	"testing"
)

func intersectionStringArrays(A []string, B []string) []string {
	intersection := make([]string, 0)
	for i := 0; i < len(A); i++ {
		for j := 0; j < len(B); j++ {
			test := A[i] == B[j]
			if test {
				intersection = append(intersection, A[i])
			}
		}
	}
	return intersection
}

func TestEncoding(t *testing.T) {

	nucs := []byte{'A', 'G', 'C', 'T', 'R', 'M', 'W', 'S', 'K', 'Y', 'V', 'H', 'D', 'B', 'N', '-', '?',
		'a', 'g', 'c', 't', 'r', 'm', 'w', 's', 'k', 'y', 'v', 'h', 'd', 'b', 'n'}

	lookupChar := make(map[byte][]string)

	lookupChar['A'] = []string{"A"}
	lookupChar['a'] = []string{"A"}
	lookupChar['C'] = []string{"C"}
	lookupChar['c'] = []string{"C"}
	lookupChar['G'] = []string{"G"}
	lookupChar['g'] = []string{"G"}
	lookupChar['T'] = []string{"T"}
	lookupChar['t'] = []string{"T"}
	lookupChar['R'] = []string{"A", "G"}
	lookupChar['r'] = []string{"A", "G"}
	lookupChar['Y'] = []string{"C", "T"}
	lookupChar['y'] = []string{"C", "T"}
	lookupChar['S'] = []string{"G", "C"}
	lookupChar['s'] = []string{"G", "C"}
	lookupChar['W'] = []string{"A", "T"}
	lookupChar['w'] = []string{"A", "T"}
	lookupChar['K'] = []string{"G", "T"}
	lookupChar['k'] = []string{"G", "T"}
	lookupChar['M'] = []string{"A", "C"}
	lookupChar['m'] = []string{"A", "C"}
	lookupChar['B'] = []string{"C", "G", "T"}
	lookupChar['b'] = []string{"C", "G", "T"}
	lookupChar['D'] = []string{"A", "G", "T"}
	lookupChar['d'] = []string{"A", "G", "T"}
	lookupChar['H'] = []string{"A", "C", "T"}
	lookupChar['h'] = []string{"A", "C", "T"}
	lookupChar['V'] = []string{"A", "C", "G"}
	lookupChar['v'] = []string{"A", "C", "G"}
	lookupChar['N'] = []string{"A", "C", "G", "T"}
	lookupChar['n'] = []string{"A", "C", "G", "T"}
	lookupChar['?'] = []string{"A", "C", "G", "T"}
	lookupChar['-'] = []string{"A", "C", "G", "T"}

	lookupByte := makeEncodingArray()

	for i := 0; i < len(nucs); i++ {
		for j := 0; j < len(nucs); j++ {
			nuc1 := nucs[i]
			nuc2 := nucs[j]

			nuc1Chars := lookupChar[nuc1]
			nuc2Chars := lookupChar[nuc2]

			byte1 := lookupByte[nuc1]
			byte2 := lookupByte[nuc2]

			byteDifferent := (byte1 & byte2) < 16
			byteSame := (byte1&8 == 8) && byte1 == byte2

			nucDifferent := len(intersectionStringArrays(nuc1Chars, nuc2Chars)) == 0
			nucSame := len(intersectionStringArrays([]string{strings.ToUpper(string(nuc1))}, []string{"A", "C", "G", "T"})) == 1 && strings.ToUpper(string(nuc1)) == strings.ToUpper(string(nuc2))

			test := byteDifferent == nucDifferent && byteSame == nucSame

			if !test {
				t.Errorf("problem in encoding test: %s %s", string(nuc1), string(nuc2))
			}
		}
	}
}

func TestDecoding(t *testing.T) {
	nucs := []byte{'A', 'G', 'C', 'T', 'R', 'M', 'W', 'S', 'K', 'Y', 'V', 'H', 'D', 'B', 'N', '-', '?',
		'a', 'g', 'c', 't', 'r', 'm', 'w', 's', 'k', 'y', 'v', 'h', 'd', 'b', 'n'}

	EA := makeEncodingArray()
	DA := makeDecodingArray()

	for _, nuc := range nucs {
		a := EA[nuc]
		b := DA[a]
		if strings.ToUpper(string(nuc)) != b {
			t.Errorf("problem in decoding test: %s", string(nuc))
		}
	}
}

func TestSNPs(t *testing.T) {
	refData := []byte(`>ref
ATGATG
`)
	queryData := []byte(
		`>Query1
ATGATG
>Query2
ATGATC
>Query3
ATTTTW
`)

	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)

	out := new(bytes.Buffer)

	err := snps(query, ref, false, false, 0.0, out)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != `query,SNPs
Query1,
Query2,G6C
Query3,G3T|A4T|G6W
` {
		t.Errorf("problem in TestSNPs()")
		fmt.Println(string(out.Bytes()))
	}
}

func TestSNPsHardGaps(t *testing.T) {
	refData := []byte(`>ref
ATGATG
`)
	queryData := []byte(
		`>Query1
--GATG
>Query2
ATGATC
>Query3
ATTTTW
`)

	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)

	out := new(bytes.Buffer)

	err := snps(query, ref, true, false, 0.0, out)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != `query,SNPs
Query1,A1-|T2-
Query2,G6C
Query3,G3T|A4T|G6W
` {
		t.Errorf("problem in TestSNPsHardGaps()")
		fmt.Println(string(out.Bytes()))
	}
}

func TestSNPsAggregate(t *testing.T) {
	refData := []byte(`>ref
ATGATG
`)
	queryData := []byte(
		`>Query1
ATGATG
>Query2
ATGATC
>Query3
ATTTTW
>Query4
ATTTTG
`)

	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)

	out := new(bytes.Buffer)

	err := snps(query, ref, false, true, 0.0, out)
	if err != nil {
		t.Error(err)
	}

	// strconv.FormatFloat(propMap[snp]/counter, 'f', 4, 64) + "\n")

	if string(out.Bytes()) != `change,proportion
G3T,0.500000000
A4T,0.500000000
G6C,0.250000000
G6W,0.250000000
` {
		t.Errorf("problem in TestSNPsAggregate()")
		fmt.Println(string(out.Bytes()))
	}
}

func TestSNPsAggregateThresh(t *testing.T) {
	refData := []byte(`>ref
ATGATG
`)
	queryData := []byte(
		`>Query1
ATGATG
>Query2
ATGATC
>Query3
ATTTTW
>Query4
ATTTTG
`)

	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)

	out := new(bytes.Buffer)

	err := snps(query, ref, false, true, 0.26, out)
	if err != nil {
		t.Error(err)
	}

	// strconv.FormatFloat(propMap[snp]/counter, 'f', 4, 64) + "\n")

	if string(out.Bytes()) != `change,proportion
G3T,0.500000000
A4T,0.500000000
` {
		t.Errorf("problem in TestSNPsAggregateThresh()")
		fmt.Println(string(out.Bytes()))
	}
}
