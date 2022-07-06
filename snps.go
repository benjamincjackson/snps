package main

import (
	"bufio"
	"errors"
	"io"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/spf13/cobra"
)

// encodedFastaRecord is a struct for one Fasta record
type encodedFastaRecord struct {
	ID          string
	Description string
	Seq         []byte
	idx         int
}

// snpLine is a struct for one Fasta record's SNPs
type snpLine struct {
	queryname string
	snps      []string
	idx       int
}

func openIn(inFile string) (*os.File, error) {
	var err error
	var f *os.File

	if inFile != "stdin" {
		f, err = os.Open(inFile)
		if err != nil {
			return f, err
		}
	} else {
		f = os.Stdin
	}

	return f, nil
}

func openOut(outFile string) (*os.File, error) {
	var err error
	var f *os.File

	if outFile != "stdout" {
		f, err = os.Create(outFile)
		if err != nil {
			return f, err
		}
	} else {
		f = os.Stdout
	}

	return f, nil
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

// makeEncodingArray returns an array whose indices are the byte representations
// of IUPAC codes and whose contents are Emmanual Paradis encodings
// Lower case nucleotides are mapped to their upper case nucleotides's encoding
func makeEncodingArrayHardGaps() []byte {
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
	byteArray['-'] = 4
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
	byteArray[4] = "-"
	byteArray[242] = "?"

	return byteArray
}

// readEncodeAlignment reads an alignment in fasta format to a channel
// of encodedFastaRecord structs - converting sequence to EP's bitwise coding scheme
func readEncodeAlignment(r io.Reader, hardGaps bool, chnl chan encodedFastaRecord, chnlerr chan error, cdone chan bool) {

	var encoding []byte
	switch hardGaps {
	case true:
		encoding = makeEncodingArrayHardGaps()
	case false:
		encoding = makeEncodingArray()
	}

	s := bufio.NewScanner(r)

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
			for i := range line {
				encodedLine[i] = encoding[line[i]]
			}
			seqBuffer = append(seqBuffer, encodedLine...)
		}
	}

	fr := encodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, idx: counter}
	chnl <- fr

	if s.Err() != nil {
		chnlerr <- s.Err()
	}

	cdone <- true
}

// getSNPs gets the SNPs between the reference and each Fasta record at a time
func getSNPs(refSeq []byte, cFR chan encodedFastaRecord, cSNPs chan snpLine, cErr chan error) {

	DA := makeDecodingArray()

	for FR := range cFR {
		SL := snpLine{}
		SL.queryname = FR.ID
		SL.idx = FR.idx
		SNPs := make([]string, 0)
		for i, nuc := range FR.Seq {
			if (refSeq[i] & nuc) < 16 {
				snpLine := DA[refSeq[i]] + strconv.Itoa(i+1) + DA[nuc]
				SNPs = append(SNPs, snpLine)
			}
		}
		SL.snps = SNPs
		cSNPs <- SL
	}

	return
}

// writeOutput writes the output to stdout as it arrives. It uses a map to write things
// in the same order as they are in the input file.
func writeOutput(w io.Writer, cSNPs chan snpLine, cErr chan error, cWriteDone chan bool) {

	outputMap := make(map[int]snpLine)

	counter := 0

	var err error

	_, err = w.Write([]byte("query,SNPs\n"))
	if err != nil {
		cErr <- err
		return
	}

	for snpLine := range cSNPs {

		outputMap[snpLine.idx] = snpLine

		for {
			if SL, ok := outputMap[counter]; ok {
				_, err = w.Write([]byte(SL.queryname + "," + strings.Join(SL.snps, "|") + "\n"))
				if err != nil {
					cErr <- err
					return
				}
				delete(outputMap, counter)
				counter++
			} else {
				break
			}
		}
	}

	cWriteDone <- true
}

func aggregateWriteOutput(w io.Writer, threshold float64, cSNPs chan snpLine, cErr chan error, cWriteDone chan bool) {

	propMap := make(map[string]float64)

	var err error

	_, err = w.Write([]byte("change,proportion\n"))
	if err != nil {
		cErr <- err
	}

	counter := 0.0

	for snpLine := range cSNPs {
		counter++
		for _, snp := range snpLine.snps {
			if _, ok := propMap[snp]; ok {
				propMap[snp]++
			} else {
				propMap[snp] = 1.0
			}
		}
	}

	order := make([]string, 0)
	for k, _ := range propMap {
		order = append(order, k)
	}

	sort.SliceStable(order, func(i, j int) bool {
		pos_i, err := strconv.Atoi(order[i][1 : len(order[i])-1])
		if err != nil {
			cErr <- err
		}
		pos_j, err := strconv.Atoi(order[j][1 : len(order[j])-1])
		if err != nil {
			cErr <- err
		}
		alt_i := order[i][len(order[i])-1]
		alt_j := order[j][len(order[j])-1]
		return pos_i < pos_j || (pos_i == pos_j && alt_i < alt_j)
	})

	for _, snp := range order {
		if propMap[snp]/counter < threshold {
			continue
		}
		_, err = w.Write([]byte(snp + "," + strconv.FormatFloat(propMap[snp]/counter, 'f', 9, 64) + "\n"))
		if err != nil {
			cErr <- err
		}
	}

	cWriteDone <- true
}

// Run the program
func snps(rQ io.Reader, rR io.Reader, hardGaps bool, aggregate bool, threshold float64, w io.Writer) error {

	cErr := make(chan error)

	cRef := make(chan encodedFastaRecord)
	cRefDone := make(chan bool)

	cFR := make(chan encodedFastaRecord)
	cFRDone := make(chan bool)

	cSNPs := make(chan snpLine, runtime.NumCPU())
	cSNPsDone := make(chan bool)

	cWriteDone := make(chan bool)

	go readEncodeAlignment(rR, hardGaps, cRef, cErr, cRefDone)

	var refSeq []byte

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case FR := <-cRef:
			refSeq = FR.Seq
		case <-cRefDone:
			close(cRef)
			n--
		}
	}

	go readEncodeAlignment(rQ, hardGaps, cFR, cErr, cFRDone)

	switch aggregate {
	case true:
		go aggregateWriteOutput(w, threshold, cSNPs, cErr, cWriteDone)
	case false:
		go writeOutput(w, cSNPs, cErr, cWriteDone)
	}

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
		cSNPsDone <- true
	}()

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cFRDone:
			close(cFR)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cSNPsDone:
			close(cSNPs)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cWriteDone:
			n--
		}
	}

	return nil
}

var snpsReference string
var snpsQuery string
var snpsOutfile string
var hardGaps bool
var aggregate bool
var thresh float64

func init() {
	mainCmd.Flags().StringVarP(&snpsReference, "reference", "r", "", "Reference sequence, in fasta format")
	mainCmd.Flags().StringVarP(&snpsQuery, "query", "q", "stdin", "Alignment of sequences to find snps in, in fasta format")
	mainCmd.Flags().StringVarP(&snpsOutfile, "outfile", "o", "stdout", "Output to write")
	mainCmd.Flags().BoolVarP(&hardGaps, "hard-gaps", "", false, "don't treat alignment gaps as missing data")
	mainCmd.Flags().BoolVarP(&aggregate, "aggregate", "", false, "report the proportions of each change")
	mainCmd.Flags().Float64VarP(&thresh, "threshold", "", 0.0, "if --aggregate, only report snps with a freq above this value")

	mainCmd.Flags().Lookup("hard-gaps").NoOptDefVal = "true"
	mainCmd.Flags().Lookup("aggregate").NoOptDefVal = "true"

	mainCmd.Flags().SortFlags = false
}

var mainCmd = &cobra.Command{
	Use:   "snps",
	Short: "snps...",
	Long:  `snps...`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		queryIn, err := openIn(snpsQuery)
		if err != nil {
			return err
		}
		defer queryIn.Close()

		refIn, err := openIn(snpsReference)
		if err != nil {
			return err
		}
		defer refIn.Close()

		snpsOut, err := openOut(snpsOutfile)
		if err != nil {
			return err
		}
		defer snpsOut.Close()

		err = snps(queryIn, refIn, hardGaps, aggregate, thresh, snpsOut)

		return err
	},
}

func main() {
	mainCmd.Execute()
}
