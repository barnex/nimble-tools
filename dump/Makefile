all: gofmt 6g #gccgo

6g:
	go install -v 

gofmt:
	gofmt -w *.go

GCCGO=gccgo -gccgoflags '-static-libgcc -O3'

gccgo:
	go build -v -compiler $(GCCGO)

