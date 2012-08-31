package main

// Output for gnuplot's "splot"
// Author: Arne Vansteenkiste

import (
	"bufio"
	"fmt"
	"nimble-cube/core"
	"nimble-cube/dump"
	"os"
)

func dumpGnuplot(f *dump.Frame, file string) {

	out, err := os.OpenFile(file, os.O_CREATE|os.O_TRUNC|os.O_WRONLY, 0666)

	out_buffered := bufio.NewWriter(out)

	core.Fatal(err)
	defer out.Close()

	data := f.Tensors()
	gridsize := f.Size()[1:]
	cellsize := f.MeshStep
	ncomp := len(data)

	// Here we loop over X,Y,Z, not Z,Y,X, because
	// internal in C-order == external in Fortran-order
	for i := 0; i < gridsize[X]; i++ {
		x := float64(i) * cellsize[X]
		for j := 0; j < gridsize[Y]; j++ {
			y := float64(j) * cellsize[Y]
			for k := 0; k < gridsize[Z]; k++ {
				z := float64(k) * cellsize[Z]
				_, err := fmt.Fprint(out, z, " ", y, " ", x, "\t")
				core.Fatal(err)
				for c := 0; c < ncomp; c++ {
					_, err := fmt.Fprint(out_buffered, data[core.SwapIndex(c, ncomp)][i][j][k], " ") // converts to user space.
					core.Fatal(err)
				}
				_, err = fmt.Fprint(out, "\n")
				core.Fatal(err)
			}
			_, err := fmt.Fprint(out, "\n")
			core.Fatal(err)
		}
		core.Fatal(err)
	}
}
