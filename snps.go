package main

import (
	"log"
	"bufio"
	"os"
	"strings"
	"errors"
	"runtime"
	"sync"
	"strconv"
)

// encodedFastaRecord is a struct for one Fasta record
type encodedFastaRecord struct {
	ID          string
	Description string
	Seq         []byte
	idx int
}

// snpLine is a struct for one Fasta record's SNPs
type snpLine struct {
	queryname string
	snps []string
	idx int
}

// makeEncodingArray returns an array whose indices are the byte representations
// of IUPAC codes and whose contents are Emmanual Paradis encodings
// Lower case nucleotides are mapped to their upper case nucleotides's encoding
func makeEncodingArray() []byte {
	byteArray := make([]byte, 256)

	byteArray['A'] = 136
	byteArray['a'] = 136
	byteArray['G'] = 72
	byteArray['g'] = 72
	byteArray['C'] = 40
	byteArray['c'] = 40
	byteArray['T'] = 24
	byteArray['t'] = 24
	byteArray['R'] = 192
	byteArray['r'] = 192
	byteArray['M'] = 160
	byteArray['m'] = 160
	byteArray['W'] = 144
	byteArray['w'] = 144
	byteArray['S'] = 96
	byteArray['s'] = 96
	byteArray['K'] = 80
	byteArray['k'] = 80
	byteArray['Y'] = 48
	byteArray['y'] = 48
	byteArray['V'] = 224
	byteArray['v'] = 224
	byteArray['H'] = 176
	byteArray['h'] = 176
	byteArray['D'] = 208
	byteArray['d'] = 208
	byteArray['B'] = 112
	byteArray['b'] = 112
	byteArray['N'] = 240
	byteArray['n'] = 240
	byteArray['-'] = 244
	byteArray['?'] = 242

	return byteArray
}

// makeDecodingArray returns an array whose indices are Emmanual Paradis encodings
// of IUPAC codes and whose contents are IUPAC codes as strings
func makeDecodingArray() []string {
	byteArray := make([]string, 256)

	byteArray[136] = "A"
	byteArray[72] = "G"
	byteArray[40] = "C"
	byteArray[24] = "T"
	byteArray[192] = "R"
	byteArray[160] = "M"
	byteArray[144] = "W"
	byteArray[96] = "S"
	byteArray[80] = "K"
	byteArray[48] = "Y"
	byteArray[224] = "V"
	byteArray[176] = "H"
	byteArray[208] = "D"
	byteArray[112] = "B"
	byteArray[240] = "N"
	byteArray[244] = "-"
	byteArray[242] = "?"

	return byteArray
}

// readEncodeAlignment reads an alignment in fasta format to a channel
// of encodedFastaRecord structs - converting sequence to EP's bitwise coding scheme
func readEncodeAlignment(infile string, chnl chan encodedFastaRecord, chnlerr chan error, cdone chan bool) {

	f, err := os.Open(infile)

	if err != nil {
		chnlerr <- err
	}

	defer f.Close()

	encoding := makeEncodingArray()
	// encoding := MakeByteDict2()

	s := bufio.NewScanner(f)

	first := true

	var id string
	var description string
	var seqBuffer []byte
	var line []byte

	counter := 0

	for s.Scan() {
		line = s.Bytes()

		if first {

			if line[0] != '>' {
				chnlerr <- errors.New("badly formatted fasta file")
			}

			description = string(line[1:])
			id = strings.Fields(description)[0]

			first = false

		} else if line[0] == '>' {

			fr := encodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, idx: counter}
			chnl <- fr
			counter++

			description = string(line[1:])
			id = strings.Fields(description)[0]
			seqBuffer = make([]byte, 0)

		} else {
			encodedLine := make([]byte, len(line))
			for i := range(line) {
				encodedLine[i] = encoding[line[i]]
			}
			seqBuffer = append(seqBuffer, encodedLine...)
		}
	}

	fr := encodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, idx: counter}
	chnl <- fr

	err = s.Err()
	if err != nil {
		chnlerr <- err
	}

	cdone <- true
}

// getSNPs gets the SNPs between the reference and each Fasta record at a time
func getSNPs(refSeq []byte, cFR chan encodedFastaRecord, cSNPs chan snpLine, cErr chan error) {

	DA := makeDecodingArray()

	for FR := range(cFR) {
		SL := snpLine{}
		SL.queryname = FR.ID
		SL.idx = FR.idx
		SNPs := make([]string, 0)
		for i, nuc := range(FR.Seq) {
			if (refSeq[i] & nuc) < 16 {
				snpLine := DA[refSeq[i]] + strconv.Itoa(i + 1) + DA[nuc]
				SNPs = append(SNPs, snpLine)
			}
		}
		SL.snps = SNPs
		cSNPs<- SL
	}

	return
}

// writeOutput writes the output to stdout as it arrives. It uses a map to write things
// in the same order as they are in the input file.
func writeOutput(cSNPs chan snpLine, cErr chan error, cWriteDone chan bool) {

	outputMap := make(map[int]snpLine)

	counter := 0

	f := os.Stdout
	defer f.Close()

	_, err := f.WriteString("query,SNPs\n")
	if err != nil {
		cErr <- err
	}

	for snpLine := range cSNPs {
		outputMap[snpLine.idx] = snpLine

		if SL, ok := outputMap[counter]; ok {
			_, err := f.WriteString(SL.queryname + "," + strings.Join(SL.snps, "|") + "\n")
			if err != nil {
				cErr <- err
			}
			delete(outputMap, counter)
			counter++
		} else {
			continue
		}
	}

	for n := 1; n > 0; {
		if len(outputMap) == 0 {
			n--
			break
		}
		SL := outputMap[counter]
		_, err := f.WriteString(SL.queryname + "," + strings.Join(SL.snps, "|") + "\n")
		if err != nil {
			cErr <- err
		}
		delete(outputMap, counter)
		counter++
	}

	cWriteDone <- true
}

// Run the program
func main() {

	if len(os.Args) != 3 {
		os.Stderr.WriteString("Usage: ./snps reference.fasta alignment.fasta > snps.csv")
	} else {
		referenceFile := os.Args[1]
		alignmentFile := os.Args[2]

		cErr := make(chan error)

		cRef := make(chan encodedFastaRecord)
		cRefDone := make(chan bool)

		cFR := make(chan encodedFastaRecord)
		cFRDone := make(chan bool)

		cSNPs := make(chan snpLine, runtime.NumCPU())
		cSNPsDone := make(chan bool)


		cWriteDone := make(chan bool)

		go readEncodeAlignment(referenceFile, cRef, cErr, cRefDone)

		var refSeq []byte

		for n := 1; n > 0; {
			select {
			case err := <-cErr:
				log.Fatal(err)
			case FR := <-cRef:
				refSeq = FR.Seq
			case <-cRefDone:
				close(cRef)
				n--
			}
		}

		go readEncodeAlignment(alignmentFile, cFR, cErr, cFRDone)

		go writeOutput(cSNPs, cErr, cWriteDone)

		var wgSNPs sync.WaitGroup
		wgSNPs.Add(runtime.NumCPU())

		for n := 0; n < runtime.NumCPU(); n++ {
			go func() {
				getSNPs(refSeq, cFR, cSNPs, cErr)
				wgSNPs.Done()
			}()
		}

		go func() {
			wgSNPs.Wait()
			cSNPsDone<- true
		}()

		for n := 1; n > 0; {
			select {
			case err := <-cErr:
				log.Fatal(err)
			case <-cFRDone:
				close(cFR)
				n--
			}
		}

		for n := 1; n > 0; {
			select {
			case err := <-cErr:
				log.Fatal(err)
			case <-cSNPsDone:
				close(cSNPs)
				n--
			}
		}

		for n := 1; n > 0; {
			select {
			case err := <-cErr:
				log.Fatal(err)
			case <-cWriteDone:
				n--
			}
		}
	}
}
