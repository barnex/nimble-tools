package main

// Calculate dispersion relations from 4D time+space FFT.
// Author: Arne Vansteenkiste.

import (
	"flag"
	"fmt"
	"github.com/barnex/fftw3"
	"github.com/barnex/fftw3/fftwf"
	"github.com/barnex/reshape"
	"math"
	"nimble-cube/core"
	"nimble-cube/dump"
	"os"
)

func main() {
	flag.Parse()
	data, logicSize := loadData()
	list := reshape.ContiguousR4(data)
	fft := fftwf.PlanDFTR2C(logicSize[:], list, reshape.RtoC(list), fftwf.ESTIMATE)
	say("windowing...")
	window(data)
	say("FFT...")
	fft.Execute()
	fSpectrum(data)
	rmDC(data)
	dispersionX(data)
	dispersionY(data)
	dispersionZ(data)
	dispersion(data)
}

func window(data [][][][]float32) {
	N := len(data)
	for n := range data {
		for x := range data[n] {
			for y := range data[0][0] {
				for z := range data[0][0][0] {
					data[n][x][y][z] *= hanning(n, N)
				}
			}
		}
	}
}

func hanning(n, N int) float32 {
	return float32(0.5 * (1 - math.Cos((2*math.Pi*float64(n))/float64(N-1))))
}

func rmDC(data [][][][]float32) {
	for x := range data[0] {
		for y := range data[0][0] {
			for z := range data[0][0][0] {
				data[0][x][y][z] = 0
			}
		}
	}
}

func fSpectrum(data [][][][]float32) {
	say("frequencyspectrum")
	out, err := os.OpenFile("frequencyspectrum", os.O_CREATE|os.O_TRUNC|os.O_WRONLY, 0666)
	core.Fatal(err)
	defer out.Close()

	for i := 1; i < len(data); i++ {
		k := data[i]
		power := 0.
		K := reshape.RtoC(reshape.ContiguousR3(k))
		for _, s := range K {
			power += float64(real(s)*real(s) + imag(s)*imag(s))
		}
		f := float64(i) / totaltime
		fmt.Fprintln(out, f, "\t", power)
	}
}

func dispersion(data [][][][]float32) {
	say("dispersion")
	out, err := os.OpenFile("dispersion", os.O_CREATE|os.O_TRUNC|os.O_WRONLY, 0666)
	core.Fatal(err)
	defer out.Close()

	nx, ny, nz := len(data[0]), len(data[0][0]), len(data[0][0][0])

	for f := 1; f < len(data)/2; f++ {
		k := data[f]
		freq := float64(f) / totaltime
		power := make([]float64, nx+ny+nz)

		for x := range k {
			for y := range k[x] {
				for z := range k[x][y] {
					power[int(math.Sqrt(float64(x*x+y*y+z*z)))] += float64(k[x][y][z] * k[x][y][z])
				}
			}
		}

		for i := range power {
			wvec := float64(i) / meshsize[0] // !!!!!!! Assumes cubic cells.
			fmt.Fprintln(out, freq, wvec, power[i])
		}
		fmt.Fprintln(out)
	}
}

type my [][2]float64

func (m *my) Len() int           { return len(*m) }
func (m *my) Less(i, j int) bool { return (*m)[i][0] < (*m)[j][0] }
func (m *my) Swap(i, j int)      { (*m)[i], (*m)[j] = (*m)[j], (*m)[i] }

func dispersionX(data [][][][]float32) {
	say("dispersionx")
	out, err := os.OpenFile("dispersionx", os.O_CREATE|os.O_TRUNC|os.O_WRONLY, 0666)
	core.Fatal(err)
	defer out.Close()

	for f := 1; f < len(data)/2; f++ {
		k := data[f]
		freq := float64(f) / totaltime

		for x := range k {
			wvec := float64(x) / meshsize[0]
			power := 0.
			for y := range k[x] {
				for z := range k[x][y] {
					power += float64(k[x][y][z] * k[x][y][z])
				}
			}
			fmt.Fprintln(out, freq, wvec, power)
		}
		fmt.Fprintln(out)
	}
}

func dispersionY(data [][][][]float32) {
	say("dispersiony")
	out, err := os.OpenFile("dispersiony", os.O_CREATE|os.O_TRUNC|os.O_WRONLY, 0666)
	core.Fatal(err)
	defer out.Close()

	for f := 1; f < len(data)/2; f++ {
		k := data[f]
		freq := float64(f) / totaltime
		for y := range k[0] {
			wvec := float64(y) / meshsize[1]
			power := 0.
			for x := range k {
				for z := range k[0][0] {
					power += float64(k[x][y][z] * k[x][y][z])
				}
			}
			fmt.Fprintln(out, freq, wvec, power)
		}
		fmt.Fprintln(out)
	}
}

func dispersionZ(data [][][][]float32) {
	say("dispersionz")
	out, err := os.OpenFile("dispersionz", os.O_CREATE|os.O_TRUNC|os.O_WRONLY, 0666)
	core.Fatal(err)
	defer out.Close()

	for f := 1; f < len(data)/2; f++ {
		k := data[f]
		freq := float64(f) / totaltime
		for z := range k[0][0] {
			wvec := float64(z) / meshsize[2]
			power := 0.
			for x := range k {
				for y := range k[0] {
					power += float64(k[x][y][z] * k[x][y][z])
				}
			}
			fmt.Fprintln(out, freq, wvec, power)
		}
		fmt.Fprintln(out)
	}
}

var totaltime = 0.
var meshsize [3]float64

func loadData() (array [][][][]float32, logicsize [4]int) {
	N := flag.NArg() // Number of dump files = Number of time points
	say("loading", N, "input files...")
	var list []float32

	frames := dump.ReadAllFiles(flag.Args(), dump.CRC_ENABLED)
	var size [3]int
	t := 0 // time
	for f := range frames {
		// init sizes from first frame.
		if size[0] == 0 {
			size = f.MeshSize
			say("mesh size is", size)
			outSize := fftw3.R2COutputSizeFloats(size[:])
			store := [4]int{N, outSize[0], outSize[1], outSize[2]}
			list = make([]float32, prod(store[:]))
			say("list size is", len(list))
			array = reshape.R4(list, store)
			logicsize = [4]int{N, f.MeshSize[0], f.MeshSize[1], f.MeshSize[2]}
			for i := range meshsize {
				meshsize[i] = f.MeshStep[i] * float64(f.MeshSize[i])
			}
		}
		totaltime = f.Time

		in := f.Tensors()[0]

		for i := range in {
			for j := range in[i] {
				for k := range in[i][j] {
					array[t][i][j][k] = in[i][j][k]
				}
			}
		}
		t++
	}
	say("data loaded")
	return
}

func say(msg ...interface{}) {
	fmt.Fprintln(os.Stderr, msg...)
}

func prod(n []int) int {
	prod := 1
	for i := range n {
		prod *= n[i]
	}
	return prod
}
