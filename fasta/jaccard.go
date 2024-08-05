package main

import (
	"bufio"
	"fmt"
	"os"

	"github.com/akamensky/argparse"
	"gonum.org/v1/gonum/stat/combin"
)

func Shingle(s string, k int) []string {

	/*
		split string into unique n-grams
	*/

	out := make([]string, 0, len(s)) //max is the number of bytes (you get with n-gram 1)
	m := make(map[string]int)

	//split string into subsequent runes - can also be bytes in this case but anyway
	runeS := []rune(s)
	//do not repeat if occurring multiple time - we need a map
	for i := 0; i < len(runeS)-k+1; i++ {
		//m[string(runeS[i:i+k])]++ //we may also want to count? unique elenments? Not needed for the time being
		k := string(runeS[i : i+k])
		_, ok := m[k]
		if !ok {
			m[k]++
			out = append(out, k)
		}
	}

	return out
}

func Union(s1, s2 []string) [][]rune {

	/*
		union of 2 shingles
	*/

	m := make(map[string]bool)

	for _, item := range s1 {
		m[item] = true
	}

	for _, item := range s2 {
		if _, ok := m[item]; !ok {
			//re-allocating?
			s1 = append(s1, item)
		}
	}

	// convert a to rune matrix (with x -> words and y -> characters)
	out := make([][]rune, len(s1))

	for i, word := range s1 {
		out[i] = []rune(word)
	}

	return out

}

func Jaccard(s1, s2 string, nmer int) float32 {

	/*
		calculate actual jaccard similarity
	*/

	//if any of the strings is empty or nmer less than 1, return 0
	if s1 == "" || s2 == "" || nmer < 1 {
		return 0
	}

	splts1 := Shingle(s1, nmer)
	splts2 := Shingle(s2, nmer)

	// convert to rune array
	rs1 := make([][]rune, len(splts1))
	for i, str := range splts1 {
		rs1[i] = []rune(str)
	}
	rs2 := make([][]rune, len(splts2))
	for i, str := range splts2 {
		rs2[i] = []rune(str)
	}

	// make union and get length
	un := len(Union(splts1, splts2))

	//jaccard
	//length of intersection divided by length of union
	in := len(rs1) + len(rs2) - un

	return float32(in) / float32(un)
}

func RFasta(fafile string) map[string]string {

	set := make(map[string]string)
	fa, _ := os.Open(fafile)
	defer fa.Close()

	b := bufio.NewScanner(fa)
	s := ""
	b.Scan()
	line := b.Text()
	id := line[1:]
	set[id] = ""

	for b.Scan() {

		line := b.Text()

		if line[0] == '>' {

			s = ""
			id = line[1:]
			set[id] = s

		} else {

			set[id] += line

		}

	}

	return set

}

func main() {

	//small use case
	//s1 := "abcde"
	//s2 := "abdcde"
	//ngram := 2

	/*
		//test shingle
		ss1 := Shingle(s1, ngram)
		ss2 := Shingle(s2, ngram)
		fmt.Println(ss1, ss2)

		//test union
		fmt.Println(Union(ss1, ss2))

		//test jaccard
		jacc := Jaccard(s1, s2, ngram)
		fmt.Println(jacc)
	*/

	parser := argparse.NewParser("jaccard", "calculate jaccard index")

	f := parser.String("f", "fasta", &argparse.Options{Required: true, Help: "Sequences in .fasta"})
	n := parser.Int("n", "ngram", &argparse.Options{Required: false, Help: "N-gram size for jaccard index calculation", Default: 2})

	err := parser.Parse(os.Args)

	if err != nil {

		fmt.Print(parser.Usage(err))
		os.Exit(1)

	}

	famap := RFasta(*f) //read fasta
	sli := make([]string, 0, len(famap))

	for k := range famap {

		sli = append(sli, k)

	}

	num := len(sli)              //number of words
	k := 2                       //pairwise distances
	hom := make(map[string]bool) //one with respect to itseld

	gen := combin.NewCombinationGenerator(num, k)
	fmt.Print("s1\ts2\tjaccard_index\n")

	for gen.Next() {

		i1, i2 := gen.Combination(nil)[0], gen.Combination(nil)[1]
		s1 := famap[sli[i1]]
		s2 := famap[sli[i2]]

		/*
			ss1 := Shingle(s1, ngram)
			ss2 := Shingle(s2, ngram)
			fmt.Println(ss1, ss2)

			fmt.Println(Union(ss1, ss2))

			jacc := Jaccard(s1, s2, ngram)
			fmt.Println(jacc)
		*/

		_, ok := hom[s1]

		if !ok {

			hom[s1] = true
			fmt.Print(sli[i1], "\t", sli[i1], "\t", Jaccard(s1, s1, *n), "\n")

		}

		_, ok = hom[s2]

		if !ok {

			hom[s2] = true
			fmt.Print(sli[i2], "\t", sli[i2], "\t", Jaccard(s2, s2, *n), "\n")

		}

		fmt.Print(sli[i1], "\t", sli[i2], "\t", Jaccard(s1, s2, *n), "\n")

	}
}
